#include "metadata.hpp"
#include "ldzipmatrix.hpp"
#include "ldzipcompressor.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <sstream>

namespace ldzip {

template <typename T>
void copy_binary_file_with_offset(const std::string& in_file,
                                  const std::string& out_file,
                                  std::streamoff read_pos,
                                  T add_value,
                                  size_t byte_size,
                                  size_t n,
                                  size_t chunk_size = 1000000) {
    std::ifstream in(in_file, std::ios::binary);
    if (!in) throw std::runtime_error("Failed to open input file: " + in_file);
    std::cout<<" copying "<<in_file<<" to "<<out_file<<std::endl;
    if (read_pos > 0) {
        in.seekg(read_pos, std::ios::beg);
    }

    std::ofstream out(out_file, std::ios::binary | std::ios::app);
    if (!out) throw std::runtime_error("Failed to open output file: " + out_file);

    size_t remaining = n;
    if (n == 0) {
        std::ifstream size_in(in_file, std::ios::binary | std::ios::ate);
        std::streamoff bytes = size_in.tellg() - read_pos;
        if (bytes < 0 || bytes % byte_size != 0)
            throw std::runtime_error("Invalid file size for element type");

        remaining = bytes / byte_size;
    }

    while (remaining > 0) {
        size_t to_read = (chunk_size < remaining) ? chunk_size : remaining;
        std::vector<T> buffer(to_read);

        in.read(reinterpret_cast<char*>(buffer.data()), to_read * byte_size);

        if constexpr (std::is_same_v<T, ldzip::COO::Entry>) {
            for (auto& e : buffer) {
                e.j += add_value.j;
            }
        } else {
            if (add_value != T{}) {
                for (auto& v : buffer) {
                    v += add_value;
                }
            }
        }

        out.write(reinterpret_cast<char*>(buffer.data()), buffer.size() * byte_size);
        remaining -= to_read;
    }
}


void merge_vars_files(const std::vector<std::string> &prefixes, const std::string &out_file) {
    std::remove(out_file.c_str());
    std::ofstream out(out_file);
    if (!out) throw std::runtime_error("Cannot open output file: " + out_file);

    bool firstFile = true;

    for (const auto &prefix : prefixes) {
        std::ifstream in(prefix + ".vars.txt");
        if (!in) throw std::runtime_error("Cannot open input file: " + prefix + ".vars.txt");

        std::string line;

        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '#')
            {
                if (firstFile) out << line << '\n';
                continue;
            }
            out << line << '\n';
        }
        firstFile = false;
    }
}


// ---------------------------------------------------------------------------
// Helpers shared by overlapping and non-overlapping concat
// ---------------------------------------------------------------------------

// Read non-comment lines from a vars file (each line = one variant).
// Used only by merge_vars_files_overlapping to reproduce the original text.
static std::vector<std::string> read_var_lines(const std::string &vars_file) {
    std::ifstream in(vars_file);
    if (!in) throw std::runtime_error("Cannot open vars file: " + vars_file);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '#')
            lines.push_back(line);
    }
    return lines;
}

// Load Variant objects for a chunk via LDZipMatrix::readVariants.
static std::vector<Variant> load_chunk_variants(const std::string &prefix) {
    LDZipMatrix m(prefix);
    m.readVariants(prefix + ".vars.txt");
    return m.variants();
}

// Count how many trailing variants of v1 match the leading variants of v2 by ID.
static uint32_t compute_overlap_size(const std::vector<Variant> &v1,
                                     const std::vector<Variant> &v2) {
    uint32_t max_k = static_cast<uint32_t>(std::min(v1.size(), v2.size()));
    for (uint32_t k = max_k; k > 0; --k) {
        bool match = true;
        for (uint32_t i = 0; i < k && match; ++i)
            if (v1[v1.size() - k + i].id() != v2[i].id())
                match = false;
        if (match) return k;
    }
    return 0;
}

// Per-chunk [start, end) range in global merged coordinates.
// Consecutive chunks overlap when bounds[i].end > bounds[i+1].start.
struct ChunkBounds {
    uint32_t start; // first variant's global index
    uint32_t end;   // one past the last variant's global index
};

// Read all chunk variant lists, compute pairwise overlap sizes, and return
// per-chunk ChunkBounds where end[i] > start[i+1] for each overlapping pair.
static std::vector<ChunkBounds> get_overlapping_variants(
        const std::vector<std::string>         &prefixes,
        std::vector<std::vector<Variant>>      &out_variants,
        std::vector<uint32_t>                  &out_overlaps) {
    size_t n = prefixes.size();
    out_variants.resize(n);
    for (size_t ci = 0; ci < n; ++ci)
        out_variants[ci] = load_chunk_variants(prefixes[ci]);

    out_overlaps.assign(n > 1 ? n - 1 : 0, 0);
    for (size_t ci = 0; ci + 1 < n; ++ci) {
        out_overlaps[ci] = compute_overlap_size(out_variants[ci], out_variants[ci + 1]);
        std::cout << " Overlap between chunk " << ci << " and chunk " << (ci + 1)
                  << ": " << out_overlaps[ci] << " variants\n";
    }

    std::vector<ChunkBounds> bounds(n);
    bounds[0] = {0, static_cast<uint32_t>(out_variants[0].size())};
    for (size_t ci = 1; ci < n; ++ci) {
        bounds[ci].start = bounds[ci - 1].end - out_overlaps[ci - 1];
        bounds[ci].end   = bounds[ci].start + static_cast<uint32_t>(out_variants[ci].size());
    }
    return bounds;
}

// Write merged vars file, skipping the leading overlap variants of each
// non-first chunk (those were already written from the previous chunk).
static void merge_vars_files_overlapping(
        const std::vector<std::string> &prefixes,
        const std::vector<std::vector<std::string>> &all_var_lines,
        const std::vector<uint32_t> &overlap_sizes,
        const std::string &out_file) {
    std::remove(out_file.c_str());
    std::ofstream out(out_file);
    if (!out) throw std::runtime_error("Cannot open output vars file: " + out_file);

    // Copy header from the first chunk
    {
        std::ifstream in(prefixes[0] + ".vars.txt");
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '#') { out << line << '\n'; break; }
        }
    }

    size_t n_chunks = prefixes.size();
    for (size_t ci = 0; ci < n_chunks; ++ci) {
        uint32_t skip = (ci > 0) ? overlap_sizes[ci - 1] : 0;
        for (size_t j = skip; j < all_var_lines[ci].size(); ++j)
            out << all_var_lines[ci][j] << '\n';
    }
}

// ---------------------------------------------------------------------------
// Per-column helpers used by concat_ldzip_overlapping
// ---------------------------------------------------------------------------

// Validate that all chunks are mutually compatible and use FULL format.
// Overlapping concat is restricted to FULL format (tabularPlink output).
static void validate_chunks_for_overlap(const std::vector<LDZipMatrix> &inputs) {
    if (inputs[0].format() != MatrixFormat::FULL)
        throw std::runtime_error(
            "Overlapping concat requires FULL matrix format "
            "(produced by plinkTabular compression).");

    for (size_t ci = 1; ci < inputs.size(); ++ci) {
        if (inputs[ci].format() != inputs[0].format())
            throw std::runtime_error("Chunk " + std::to_string(ci) +
                                     " has a different matrix format than chunk 0.");
        if (inputs[ci].bitsEnum() != inputs[0].bitsEnum())
            throw std::runtime_error("Chunk " + std::to_string(ci) +
                                     " has a different bit width than chunk 0.");
        auto s0 = inputs[0].stats_available();
        auto si = inputs[ci].stats_available();
        if (si.size() != s0.size() || si != s0)
            throw std::runtime_error("Chunk " + std::to_string(ci) +
                                     " has different statistics than chunk 0.");
    }
}

// Translate local uint32 row indices to global size_t indices by adding an offset.
static std::vector<size_t> to_global_rows(const std::vector<uint32_t> &local_rows,
                                          uint32_t global_offset) {
    std::vector<size_t> out(local_rows.size());
    for (size_t k = 0; k < local_rows.size(); ++k)
        out[k] = local_rows[k] + global_offset;
    return out;
}

// From a decoded column, keep only entries with local row < threshold.
// Translate kept rows to global indices by adding global_offset.
static std::pair<std::vector<size_t>, std::vector<float>>
filter_above_diagonal(const std::vector<uint32_t> &rows,
                      const std::vector<float>    &xs,
                      uint32_t threshold,
                      uint32_t global_offset) {
    std::vector<size_t> out_rows;
    std::vector<float>  out_xs;
    for (size_t k = 0; k < rows.size(); ++k) {
        if (rows[k] < threshold) {
            out_rows.push_back(rows[k] + global_offset);
            out_xs.push_back(xs[k]);
        }
    }
    return {out_rows, out_xs};
}

// From a decoded column, keep only entries with local row >= threshold
// (diagonal + below-diagonal).  Translate rows to global indices.
static std::pair<std::vector<size_t>, std::vector<float>>
filter_below_diagonal(const std::vector<uint32_t> &rows,
                      const std::vector<float>    &xs,
                      uint32_t threshold,
                      uint32_t global_offset) {
    std::vector<size_t> out_rows;
    std::vector<float>  out_xs;
    for (size_t k = 0; k < rows.size(); ++k) {
        if (rows[k] >= threshold) {
            out_rows.push_back(rows[k] + global_offset);
            out_xs.push_back(xs[k]);
        }
    }
    return {out_rows, out_xs};
}

// Append `src` onto `dst` in-place.
static void append_column(std::pair<std::vector<size_t>, std::vector<float>> &dst,
                          const std::pair<std::vector<size_t>, std::vector<float>> &src) {
    dst.first.insert(dst.first.end(), src.first.begin(), src.first.end());
    dst.second.insert(dst.second.end(), src.second.begin(), src.second.end());
}

// Decode a non-overlap column from a chunk, translate rows to global coordinates,
// and write it via the compressor (ColumnStream mode, active_col must match cidx).
static void push_nonoverlap_col(LDZipCompressor       &compressor,
                                const LDZipMatrix     &chunk,
                                uint32_t               local_col,
                                uint32_t               global_offset,
                                Stat                   stat,
                                uint32_t               active_col) {
    auto local_rows  = chunk.get_i(local_col);
    auto xs          = chunk.get_x(local_col, stat);
    auto global_rows = to_global_rows(local_rows, global_offset);
    compressor.push_column_raw(active_col, xs, global_rows, stat);
}

// Decode an overlap column from two chunks, apply the diagonal split rule, merge,
// and write via the compressor.
//   chunk1 @ j1 : keep rows < j1  (above-diagonal in chunk1's local coordinates)
//   chunk2 @ j2 : keep rows >= j2 (diagonal + below-diagonal in chunk2's local coords)
// Rows from chunk1 are translated by global_off1; rows from chunk2 by global_off2.
static void push_overlap_col(LDZipCompressor   &compressor,
                             const LDZipMatrix &chunk1, const LDZipMatrix &chunk2,
                             uint32_t           j1,    uint32_t           j2,
                             uint32_t           global_off1, uint32_t global_off2,
                             Stat stat, uint32_t active_col) {
    auto rows1 = chunk1.get_i(j1);
    auto xs1   = chunk1.get_x(j1, stat);
    auto rows2 = chunk2.get_i(j2);
    auto xs2   = chunk2.get_x(j2, stat);

    auto merged = filter_above_diagonal(rows1, xs1, j1, global_off1);
    append_column(merged, filter_below_diagonal(rows2, xs2, j2, global_off2));

    compressor.push_column_raw(active_col, merged.second, merged.first, stat);
}

// ---------------------------------------------------------------------------
// concat_ldzip_overlapping
// ---------------------------------------------------------------------------
// Concatenates multiple .ldzip chunks that may have overlapping variant ranges.
// Overlap is auto-detected by comparing variant IDs in consecutive vars files.
// For each column in the overlapping region the above-diagonal entries
// (row < column) are taken from the earlier chunk and the diagonal + below-diagonal
// entries (row >= column) are taken from the later chunk.
// Restricted to FULL matrix format (produced by plinkTabular compression).
void concat_ldzip_overlapping(const std::vector<std::string> &prefixes,
                              const std::string &out_prefix) {

    if (prefixes.empty())
        throw std::runtime_error("No input chunks provided.");

    // -----------------------------------------------------------------------
    // Step 1 – read variants and detect pairwise overlaps
    // -----------------------------------------------------------------------
    std::vector<std::vector<Variant>> all_variants;
    std::vector<uint32_t>             overlap_sizes;
    auto bounds  = get_overlapping_variants(prefixes, all_variants, overlap_sizes);
    uint32_t total_n = bounds.back().end;

    // -----------------------------------------------------------------------
    // Step 2 – open inputs and validate
    // -----------------------------------------------------------------------
    size_t n_chunks = prefixes.size();
    std::vector<LDZipMatrix> inputs;
    inputs.reserve(n_chunks);
    for (auto &p : prefixes) inputs.emplace_back(p);

    validate_chunks_for_overlap(inputs);

    auto stats = inputs[0].stats_available();
    Stat stat  = stats[0];

    // -----------------------------------------------------------------------
    // Step 3 – create output compressor (ColumnStream mode)
    // -----------------------------------------------------------------------
    LDZipCompressor compressor(total_n, total_n,
                               inputs[0].format(), stats, inputs[0].bitsEnum(),
                               out_prefix, LDZipCompressor::Mode::ColumnStream);

    // -----------------------------------------------------------------------
    // Step 4 – push columns in global order
    // -----------------------------------------------------------------------
    uint32_t active_col = 0;

    for (size_t ci = 0; ci < n_chunks; ++ci) {
        uint32_t n_cols  = static_cast<uint32_t>(all_variants[ci].size());
        uint32_t ov_prev = (ci > 0)            ? overlap_sizes[ci - 1] : 0;
        uint32_t ov_next = (ci + 1 < n_chunks) ? overlap_sizes[ci]     : 0;

        // Non-overlap columns unique to this chunk
        for (uint32_t j = ov_prev; j < n_cols - ov_next; ++j)
            push_nonoverlap_col(compressor, inputs[ci], j,
                                bounds[ci].start, stat, active_col++);

        // Overlap columns shared with the next chunk
        if (ci + 1 < n_chunks) {
            uint32_t ci_next = ci + 1;
            for (uint32_t ov = 0; ov < ov_next; ++ov) {
                uint32_t j1 = n_cols - ov_next + ov; // local col in chunk ci
                uint32_t j2 = ov;                     // local col in chunk ci+1
                push_overlap_col(compressor,
                                 inputs[ci],      inputs[ci_next],
                                 j1, j2,
                                 bounds[ci].start, bounds[ci_next].start,
                                 stat, active_col++);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Step 5 – finalize output files and write vars
    // -----------------------------------------------------------------------
    compressor.stream_close();

    // Build raw text lines for the vars file (needed to reproduce original format)
    std::vector<std::vector<std::string>> all_var_lines(n_chunks);
    for (size_t ci = 0; ci < n_chunks; ++ci)
        all_var_lines[ci] = read_var_lines(prefixes[ci] + ".vars.txt");

    merge_vars_files_overlapping(prefixes, all_var_lines, overlap_sizes,
                                 out_prefix + ".vars.txt");
}

// ---------------------------------------------------------------------------

void concat_ldzip(const std::vector<std::string> &prefixes,
                        const std::string &out_prefix) {

    uint32_t row_offset = 0, total_rows = 0;
    uint64_t nnz_offset = 0, total_nnz = 0;
    Bits bits;

    for (const auto &pref : prefixes) {
        LDZipMatrix in(pref);
        total_rows += static_cast<uint32_t>(in.nrows());
        total_nnz += static_cast<uint64_t>(in.nnz());
    }

    LDZipMatrix in(prefixes[0]);
    LDZipCompressor concator(total_rows, total_rows, in.format(), in.stats_available(), in.bitsEnum(), out_prefix, LDZipCompressor::Mode::ColumnStream);
    concator.m().set_nnz(total_nnz);

    // Initialize p_file with 0
    {
        std::ofstream con_p(concator.m().pFile(), std::ios::binary | std::ios::app);
        uint64_t zero = 0;
        con_p.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
    }

    for (const auto &pref : prefixes) {
        std::cout << " reading : " << pref << "\n";
        LDZipMatrix in(pref);
        size_t nrows = in.nrows();
        bits = in.bitsEnum();
        uint64_t nnz = in.nnz();

        copy_binary_file_with_offset<uint64_t>(in.pFile(), concator.m().pFile(), sizeof(uint64_t), nnz_offset, sizeof(uint64_t), nrows);
        copy_binary_file_with_offset<int16_t>(in.iFile(), concator.m().iFile(), 0, 0, sizeof(int16_t), nnz);
        copy_binary_file_with_offset<ldzip::COO::Entry>(in.IFile(), concator.m().IFile(), 0, ldzip::COO::Entry{0, row_offset, 0}, sizeof(ldzip::COO::Entry), 0);
        for (Stat s : All_Stats()) 
            if (in.has_stat(s)) {
                if(bits!=Bits::B99)
                    copy_binary_file_with_offset<int32_t>(in.xFile(s), concator.m().xFile(s), 0L, 0, bits_to_int(bits) / 8, nnz);
                else
                    copy_binary_file_with_offset<float>(in.xFile(s), concator.m().xFile(s), 0L, 0.0f, sizeof(float), nnz);
            }
        row_offset += nrows;
        nnz_offset += nnz;
    }

    write_metadata_json(concator.m().metaFile(), concator.m().metaInfo());
    merge_vars_files(prefixes, out_prefix + ".vars.txt");
}


}