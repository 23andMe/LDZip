#include "ldzipconcatenator.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

namespace ldzip {

LDZipConcatenator::LDZipConcatenator(size_t nrows,
                                     size_t ncols,
                                     MatrixFormat format,
                                     const std::vector<Stat>& stats,
                                     Bits bits,
                                     const std::string& prefix)
    : overlap_merger_(nrows, ncols, format, stats, bits, prefix, LDZipCompressor::Mode::ColumnStream) {

    // Set up large buffers for all streams for better performance
    constexpr size_t buffer_size = 8 * 1024 * 1024;  // 8 MB buffer
    overlap_merger_.p_stream_.rdbuf()->pubsetbuf(nullptr, buffer_size);
    overlap_merger_.i_stream_.rdbuf()->pubsetbuf(nullptr, buffer_size);
    for (Stat s : All_Stats())
        if (overlap_merger_.m_.has_stat_[s])
            overlap_merger_.x_streams_[s].rdbuf()->pubsetbuf(nullptr, buffer_size);

    // Initialize p_file with 0
    write_initial_p_zero();
}

void LDZipConcatenator::write_initial_p_zero() {
    uint64_t zero = 0;
    overlap_merger_.p_stream_.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
    overlap_merger_.p_stream_.flush();  // Ensure it's written
}

void LDZipConcatenator::init_p_for_column(uint32_t col, uint64_t value) {
    overlap_merger_.m_.p_[col] = value;
}

uint64_t LDZipConcatenator::get_p_value(uint32_t col) const {
    return overlap_merger_.m_.p_[col];
}

void LDZipConcatenator::close() {
    for (Stat s : overlap_merger_.m_.stats_available_)
        overlap_merger_.x_streams_[s].close();
    overlap_merger_.i_stream_.close();
    overlap_merger_.p_stream_.close();

    // Set total nnz from tracked value
    overlap_merger_.m_.set_nnz(current_nnz_);

    write_metadata_json(overlap_merger_.m_.metaFile(), overlap_merger_.m_.metaInfo());
}

void LDZipConcatenator::push_decoded_column(uint32_t cidx) {
    EnumArray<float, Stat> Statvalues(-999.0f);
    for (size_t k = 0; k < merged_rows_.size(); ++k) {
        for (Stat s : All_Stats())
            if (overlap_merger_.m_.has_stat_[s])
                Statvalues[s] = merged_x_[s][k];

        overlap_merger_.push_value_(merged_rows_[k], cidx, Statvalues);
    }
    overlap_merger_.active_column_ = cidx;
    writeMergedColumn();
}

void LDZipConcatenator::writeMergedColumn() {
    // Update p_vector before writing out column
    overlap_merger_.m_.p_[overlap_merger_.active_column_ + 1] =
        overlap_merger_.m_.p_[overlap_merger_.active_column_] +
        overlap_merger_.m_.i_[overlap_merger_.active_column_].size();

    // Write p-value immediately to stream (already in absolute coordinates)
    uint64_t p_val = overlap_merger_.m_.p_[overlap_merger_.active_column_ + 1];
    overlap_merger_.p_stream_.write(reinterpret_cast<const char*>(&p_val), sizeof(p_val));

    // Write i and x data
    for (Stat s : All_Stats()) if (overlap_merger_.m_.has_stat_[s]) overlap_merger_.write_x(s);
    overlap_merger_.write_i();
}

void LDZipConcatenator::clear_merge_buffers(const LDZipMatrix& mat) {
    above_rows_.clear(); above_idx_.clear();
    below_rows_.clear(); below_idx_.clear();
    merged_rows_.clear();

    for (Stat s : All_Stats()) {
        if (mat.has_stat(s)) {
            above_x_[s].clear();
            below_x_[s].clear();
            merged_x_[s].clear();
        }
    }
}

void LDZipConcatenator::extract_above_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t threshold, uint32_t global_offset) {
    auto rows = mat.get_i(local_col);
    for (uint32_t k = 0; k < rows.size(); k++) {
        if (rows[k] <= threshold) {
            above_rows_.push_back(rows[k] + global_offset);
            above_idx_.push_back(k);
        }
    }

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            auto x = mat.get_x(local_col, s);
            for (uint32_t k : above_idx_) above_x_[s].push_back(x[k]);
        }
}

void LDZipConcatenator::extract_below_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t threshold, uint32_t global_offset) {
    auto rows = mat.get_i(local_col);
    for (uint32_t k = 0; k < rows.size(); k++) {
        if (rows[k] > threshold) {
            below_rows_.push_back(rows[k] + global_offset);
            below_idx_.push_back(k);
        }
    }

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            auto x = mat.get_x(local_col, s);
            for (uint32_t k : below_idx_) below_x_[s].push_back(x[k]);
        }
}

void LDZipConcatenator::merge_above_below_buffers(const LDZipMatrix& mat) {
    merged_rows_.insert(merged_rows_.end(), above_rows_.begin(), above_rows_.end());
    merged_rows_.insert(merged_rows_.end(), below_rows_.begin(), below_rows_.end());

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            merged_x_[s].insert(merged_x_[s].end(), above_x_[s].begin(), above_x_[s].end());
            merged_x_[s].insert(merged_x_[s].end(), below_x_[s].begin(), below_x_[s].end());
        }
}

void LDZipConcatenator::process_overlap_columns(const LDZipMatrix& current_mat,
                                                 const LDZipMatrix& next_mat,
                                                 const ChunkBoundary& current_chunk,
                                                 const ChunkBoundary& next_chunk) {
    uint32_t gs_current = current_chunk.global_start;
    uint32_t gs_next = next_chunk.global_start;

    // Reopen COO stream in append mode for overlap columns
    overlap_merger_.m().reopen_I_append();

    // Initialize p_[] for the first overlapping column
    uint32_t first_overlap_col = gs_current + current_chunk.second_overlap_start;
    init_p_for_column(first_overlap_col, current_nnz_);

    for (uint32_t lc_current = current_chunk.second_overlap_start; lc_current < current_chunk.n_variants; lc_current++) {
        uint32_t global_col = gs_current + lc_current;
        uint32_t lc_next    = lc_current - current_chunk.second_overlap_start;

        clear_merge_buffers(current_mat);

        // Above/diagonal from current chunk (row <= lc_current)
        extract_above_column(current_mat, lc_current, lc_current, gs_current);

        // Below diagonal from next chunk (row > lc_next)
        extract_below_column(next_mat, lc_next, lc_next, gs_next);

        // Merge above + below
        merge_above_below_buffers(current_mat);

        push_decoded_column(global_col);
    }

    // Calculate overlapping columns nnz and track it
    uint32_t last_overlap_col = gs_current + current_chunk.n_variants;
    uint64_t overlap_nnz = get_p_value(last_overlap_col) - get_p_value(first_overlap_col);
    current_nnz_ += overlap_nnz;
    current_I_offset_ = overlap_merger_.m_.I_.get_index();
}

void LDZipConcatenator::process_exclusive_columns(const LDZipMatrix& current_mat,
                                                   const ChunkBoundary& current_chunk) {
    size_t lc_start = current_chunk.first_overlap_end;
    size_t lc_end   = current_chunk.second_overlap_start;

    if (lc_start >= lc_end) return;  // No exclusive columns

    // Close COO stream before binary copy to avoid multiple file handles
    overlap_merger_.m().close_I();

    uint64_t p_start = current_mat.get_p(lc_start);
    uint64_t p_end   = current_mat.get_p(lc_end);
    size_t excl_nrows = lc_end - lc_start;
    uint64_t excl_nnz = p_end - p_start;
    auto [I_start, I_count] = get_I_range(current_mat.IIndexFile(), lc_start, lc_end);

    Bits bits = current_mat.bitsEnum();

    // Adjust offset: subtract p[lc_start] since source p-values are cumulative
    uint64_t p_offset_adjustment = current_nnz_ - p_start;
    uint64_t I_offset_adjustment = current_I_offset_ - (I_start / sizeof(uint32_t));

    // Use stream-based copying for p, i, and x files (keeps streams open)
    copy_binary_file_with_offset<uint64_t>(current_mat.pFile(), overlap_merger_.p_stream_,
        (lc_start + 1) * sizeof(uint64_t), p_offset_adjustment, sizeof(uint64_t), excl_nrows);

    // Use file-based copying for I.bin (COO overflow file - no stream for this)
    copy_binary_file_with_offset_file<uint32_t>(current_mat.IFile(), overlap_merger_.m_.IFile(),
        I_start, uint32_t{0}, sizeof(uint32_t), I_count);

    // Copy and adjust index entries (skip first boundary - already written by constructor or overlap)
    copy_binary_file_with_offset_file<uint64_t>(current_mat.IIndexFile(), overlap_merger_.m_.IIndexFile(),
        (lc_start + 1) * sizeof(uint64_t), I_offset_adjustment, sizeof(uint64_t), excl_nrows);

    copy_binary_file_with_offset<int16_t>(current_mat.iFile(), overlap_merger_.i_stream_,
        p_start * sizeof(int16_t), static_cast<int16_t>(0), sizeof(int16_t), excl_nnz);

    for (Stat s : All_Stats())
        if (current_mat.has_stat(s)) {
            if (bits != Bits::B99)
                copy_binary_file_with_offset<int32_t>(current_mat.xFile(s), overlap_merger_.x_streams_[s],
                    p_start * (bits_to_int(bits) / 8), 0, bits_to_int(bits) / 8, excl_nnz);
            else
                copy_binary_file_with_offset<float>(current_mat.xFile(s), overlap_merger_.x_streams_[s],
                    p_start * sizeof(float), 0.0f, sizeof(float), excl_nnz);
        }

    // Track nnz for exclusive columns
    current_nnz_ += excl_nnz;
    current_I_offset_ += I_count;
}

// Stream-based version - uses existing open stream
template <typename T>
void LDZipConcatenator::copy_binary_file_with_offset(const std::string& in_file,
                                                      std::fstream& out_stream,
                                                      std::streamoff read_pos,
                                                      T add_value,
                                                      size_t byte_size,
                                                      size_t n,
                                                      size_t chunk_size) {
    std::ifstream in(in_file, std::ios::binary);
    if (!in) throw std::runtime_error("Failed to open input file: " + in_file);
    std::cout << " copying " << in_file << " to stream" << std::endl;
    if (read_pos > 0) {
        in.seekg(read_pos, std::ios::beg);
    }

    size_t remaining = n;
    if (n == 0) {
        std::ifstream size_in(in_file, std::ios::binary | std::ios::ate);
        std::streamoff bytes = size_in.tellg() - read_pos;
        if (bytes < 0 || bytes % static_cast<std::streamoff>(byte_size) != 0)
            throw std::runtime_error("Invalid file size for element type");

        remaining = static_cast<size_t>(bytes) / byte_size;
    }

    while (remaining > 0) {
        size_t to_read = (chunk_size < remaining) ? chunk_size : remaining;
        std::vector<T> buffer(to_read);

        in.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(to_read * byte_size));

        if (add_value != T{}) {
            for (auto& v : buffer) {
                v += add_value;
            }
        }

        out_stream.write(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(buffer.size() * byte_size));
        remaining -= to_read;
    }
}

// File-based version - opens new file handle (for I.bin overflow file)
template <typename T>
void LDZipConcatenator::copy_binary_file_with_offset_file(const std::string& in_file,
                                                           const std::string& out_file,
                                                           std::streamoff read_pos,
                                                           T add_value,
                                                           size_t byte_size,
                                                           size_t n,
                                                           size_t chunk_size) {
    std::fstream out_stream(out_file, std::ios::out | std::ios::binary | std::ios::app);
    if (!out_stream) throw std::runtime_error("Failed to open output file: " + out_file);

    copy_binary_file_with_offset<T>(in_file, out_stream, read_pos, add_value, byte_size, n, chunk_size);
}

std::pair<uint64_t, uint64_t> LDZipConcatenator::get_I_range(const std::string& index_file, size_t a, size_t b) {
    std::ifstream index_in(index_file, std::ios::binary);
    if (!index_in) throw std::runtime_error("Cannot open I index file: " + index_file);

    uint64_t start_pos, end_pos;
    index_in.seekg(a * sizeof(uint64_t), std::ios::beg);
    index_in.read(reinterpret_cast<char*>(&start_pos), sizeof(uint64_t));
    index_in.seekg(b * sizeof(uint64_t), std::ios::beg);
    index_in.read(reinterpret_cast<char*>(&end_pos), sizeof(uint64_t));

    return {start_pos * sizeof(uint32_t), end_pos - start_pos};
}

// Explicit template instantiations for stream-based version
template void LDZipConcatenator::copy_binary_file_with_offset<uint64_t>(const std::string&, std::fstream&, std::streamoff, uint64_t, size_t, size_t, size_t);
template void LDZipConcatenator::copy_binary_file_with_offset<int16_t>(const std::string&, std::fstream&, std::streamoff, int16_t, size_t, size_t, size_t);
template void LDZipConcatenator::copy_binary_file_with_offset<int32_t>(const std::string&, std::fstream&, std::streamoff, int32_t, size_t, size_t, size_t);
template void LDZipConcatenator::copy_binary_file_with_offset<float>(const std::string&, std::fstream&, std::streamoff, float, size_t, size_t, size_t);

// Explicit template instantiations for file-based version (for I.bin)
template void LDZipConcatenator::copy_binary_file_with_offset_file<uint32_t>(const std::string&, const std::string&, std::streamoff, uint32_t, size_t, size_t, size_t);

namespace {

// Reads the next non-header, non-empty line and extracts the variant ID
// (3rd tab-separated field: CHROM, POS, ID, ... or entire line if single-column format).
bool read_next_variant_id(std::ifstream& in, std::string& id) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        // Check if line has tabs (multi-column format)
        size_t first_tab = line.find('\t');
        if (first_tab != std::string::npos) {
            // Multi-column format: extract 3rd field (ID)
            std::istringstream ss(line);
            std::string chrom, pos;
            if (std::getline(ss, chrom, '\t') &&
                std::getline(ss, pos,   '\t') &&
                std::getline(ss, id,    '\t'))
                return true;
            // If we can't parse 3 fields, fall through to single-column
        }

        // Single-column format: treat entire line as ID
        id = line;
        return true;
    }
    return false;
}

} // anonymous namespace

OverlapVariantInfo read_overlapping_variant_order(const std::vector<std::string>& var_files, bool overlapping) {
    if (var_files.empty()) throw std::runtime_error("No variant files provided");

    size_t n = var_files.size();
    OverlapVariantInfo info;
    info.chunks.resize(n);

    // --- Non-overlapping mode: treat each chunk as independent ---
    if (!overlapping) {
        size_t global_pos = 0;
        for (size_t i = 0; i < n; i++) {
            std::ifstream f(var_files[i]);
            if (!f) throw std::runtime_error("Cannot open variant file: " + var_files[i]);
            std::string id;
            size_t count = 0;
            while (read_next_variant_id(f, id)) ++count;

            // first_overlap_end=0 means no overlap with previous chunk
            // second_overlap_start=count means no overlap with next chunk
            info.chunks[i] = {count, 0, count, global_pos, global_pos + count};
            global_pos += count;
        }
        info.total_variants = global_pos;
        return info;
    }

    // --- Single file: trivial case ---
    if (n == 1) {
        std::ifstream f(var_files[0]);
        if (!f) throw std::runtime_error("Cannot open variant file: " + var_files[0]);
        std::string id;
        size_t count = 0;
        while (read_next_variant_id(f, id)) ++count;
        info.chunks[0] = {count, 0, count, 0, count};
        info.total_variants = count;
        return info;
    }

    std::ifstream first(var_files[0]);
    if (!first) throw std::runtime_error("Cannot open variant file: " + var_files[0]);
    std::ifstream second(var_files[1]);
    if (!second) throw std::runtime_error("Cannot open variant file: " + var_files[1]);

    std::vector<size_t> local_counts(n, 0);
    size_t global_cursor = 0;

    info.chunks[0].global_start      = 0;
    info.chunks[0].first_overlap_end = 0;

    for (size_t i = 0; i + 1 < n; i++) {
        // Step 1: read pivot = first variant of second (chunk[i+1])
        std::string pivot;
        if (!read_next_variant_id(second, pivot))
            throw std::runtime_error("Variant file is empty: " + var_files[i + 1]);
        ++local_counts[i + 1];

        // Step 2: read first (chunk[i]) until we reach the pivot
        std::string v;
        bool found_pivot = false;
        while (read_next_variant_id(first, v)) {
            ++local_counts[i];
            if (v == pivot) {
                found_pivot = true;
                break;
            }
        }

        size_t overlap_start_local_i;
        size_t new_variants_of_i;
        size_t overlap_count;

        if (found_pivot) {
            // v == pivot; local_counts[i] now includes the pivot itself
            overlap_start_local_i = local_counts[i] - 1;  // 0-based index of pivot in chunk[i]
            new_variants_of_i     = overlap_start_local_i - info.chunks[i].first_overlap_end;

            info.chunks[i + 1].global_start = global_cursor + new_variants_of_i;

            // Step 3: overlap loop — both streams are now positioned just past the pivot;
            // read one from each and verify they match until first is exhausted
            overlap_count = 1;  // pivot already matched
            std::string v0, v1;
            while (read_next_variant_id(first, v0)) {
                ++local_counts[i];
                if (!read_next_variant_id(second, v1))
                    throw std::runtime_error(
                        "File '" + var_files[i + 1] + "' exhausted during overlap with '" +
                        var_files[i] + "'");
                ++local_counts[i + 1];
                if (v0 != v1)
                    throw std::runtime_error(
                        "Overlap mismatch between '" + var_files[i] + "' and '" +
                        var_files[i + 1] + "': '" + v0 + "' vs '" + v1 + "'");
                ++overlap_count;
            }
        } else {
            // No overlap found - reset second file stream to beginning
            overlap_count = 0;
            overlap_start_local_i = local_counts[i];
            new_variants_of_i = local_counts[i] - info.chunks[i].first_overlap_end;
            info.chunks[i + 1].global_start = global_cursor + new_variants_of_i;
            // Reset second stream since we read the pivot but didn't match it
            second.close();
            second.open(var_files[i + 1]);
            if (!second) throw std::runtime_error("Cannot re-open variant file: " + var_files[i + 1]);
            local_counts[i + 1] = 0;
        }

        std::cout << " Overlap between chunk " << i << " and " << (i + 1)
                  << ": " << overlap_count << " variants\n";

        // Finalize chunk[i]
        info.chunks[i].n_variants           = local_counts[i];
        info.chunks[i].second_overlap_start = overlap_start_local_i;
        info.chunks[i].global_end           = global_cursor + new_variants_of_i + overlap_count;

        info.chunks[i + 1].first_overlap_end = local_counts[i + 1];  // == overlap_count

        global_cursor += new_variants_of_i + overlap_count;

        // Swap: first <- second (positioned after the overlap), open next file into second
        first = std::move(second);
        if (i + 2 < n) {
            second.open(var_files[i + 2]);
            if (!second) throw std::runtime_error("Cannot open variant file: " + var_files[i + 2]);
        }
    }

    // After the loop, first points to var_files[n-1] positioned right after the last overlap.
    // Drain the remaining tail.
    {
        std::string id;
        while (read_next_variant_id(first, id)) ++local_counts[n - 1];
    }

    size_t tail = local_counts[n - 1] - info.chunks[n - 1].first_overlap_end;
    info.chunks[n - 1].n_variants           = local_counts[n - 1];
    info.chunks[n - 1].second_overlap_start = local_counts[n - 1];  // no next overlap
    info.chunks[n - 1].global_end           = global_cursor + tail;
    info.total_variants                     = info.chunks[n - 1].global_end;

    return info;
}

} // namespace ldzip
