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
#include <algorithm>
#include <cmath>
#include <limits>
#include <climits>

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
// Helpers for overlapping concat
// ---------------------------------------------------------------------------

// Read non-comment lines from a vars file (each line = one variant)
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

// Extract the variant ID (field index 2, 0-based) from a tab-separated vars line
static std::string get_var_id(const std::string &line) {
    std::istringstream iss(line);
    std::string field;
    for (int i = 0; i <= 2; ++i) {
        if (!std::getline(iss, field, '\t'))
            throw std::runtime_error("Malformed vars line (fewer than 3 fields): " + line);
    }
    return field;
}

// Return the number of trailing variants of chunk[i] that match the leading
// variants of chunk[i+1] (by variant ID).  Returns 0 if there is no overlap.
static uint32_t compute_overlap_size(const std::vector<std::string> &lines1,
                                     const std::vector<std::string> &lines2) {
    uint32_t max_k = static_cast<uint32_t>(std::min(lines1.size(), lines2.size()));
    for (uint32_t k = max_k; k > 0; --k) {
        bool match = true;
        for (uint32_t i = 0; i < k && match; ++i) {
            if (get_var_id(lines1[lines1.size() - k + i]) != get_var_id(lines2[i]))
                match = false;
        }
        if (match) return k;
    }
    return 0;
}

// Write a merged vars file that deduplicates the overlap variants.
// For each chunk[ci] (ci > 0), the first overlap_sizes[ci-1] variants are
// already present from the previous chunk and are skipped.
static void merge_vars_files_overlapping(
        const std::vector<std::string> &prefixes,
        const std::vector<std::vector<std::string>> &all_var_lines,
        const std::vector<uint32_t> &overlap_sizes,
        const std::string &out_file) {
    std::remove(out_file.c_str());
    std::ofstream out(out_file);
    if (!out) throw std::runtime_error("Cannot open output vars file: " + out_file);

    // Write the header from the first chunk
    {
        std::ifstream in(prefixes[0] + ".vars.txt");
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '#') {
                out << line << '\n';
                break;
            }
        }
    }

    size_t n_chunks = prefixes.size();
    for (size_t ci = 0; ci < n_chunks; ++ci) {
        // Skip variants that are the overlap with the previous chunk
        uint32_t skip = (ci > 0) ? overlap_sizes[ci - 1] : 0;
        for (size_t j = skip; j < all_var_lines[ci].size(); ++j)
            out << all_var_lines[ci][j] << '\n';
    }
}

// Write a quantized value buffer to an output stream.
// Mirrors write_scaled_buffer in ldzipcompressor.cpp.
template <typename T>
static void write_quantized_vals(std::ofstream &out,
                                 const std::vector<float> &vals,
                                 int64_t scale) {
    std::vector<T> buf(vals.size());
    constexpr T NA_SENTINEL = std::numeric_limits<T>::min();
    for (size_t i = 0; i < vals.size(); ++i) {
        if (std::isnan(vals[i])) {
            buf[i] = NA_SENTINEL;
        } else {
            float v = std::clamp(vals[i], -1.0f, 1.0f);
            if constexpr (std::is_same_v<T, float>)
                buf[i] = v;
            else
                buf[i] = static_cast<T>(
                    std::llround(static_cast<double>(v) * static_cast<double>(scale)));
        }
    }
    out.write(reinterpret_cast<const char *>(buf.data()),
              static_cast<std::streamsize>(buf.size() * sizeof(T)));
}

// Delta-encode global_rows relative to global_col and write to i_out / coo_out.
// Also writes quantized x values for all stats to x_outs.
// Updates current_nnz and p_vec[global_col] / p_vec[global_col+1].
static void write_column_to_files(
        uint32_t global_col,
        const std::vector<uint32_t> &global_rows,   // ascending global row indices
        const EnumArray<std::vector<float>, Stat> &row_vals,
        uint64_t &current_nnz,
        std::vector<uint64_t> &p_vec,
        std::ofstream &i_out,
        COO &coo_out,
        EnumArray<std::ofstream, Stat> &x_outs,
        Bits bits,
        const std::vector<Stat> &stats_available) {

    p_vec[global_col] = current_nnz;
    size_t nnz = global_rows.size();

    if (nnz > 0) {
        // Write x values for each stat
        int64_t scale = (1LL << (bits_to_int(bits) - 1)) - 1;
        for (Stat s : stats_available) {
            switch (bits) {
                case Bits::B8:  write_quantized_vals<int8_t> (x_outs[s], row_vals[s], scale); break;
                case Bits::B16: write_quantized_vals<int16_t>(x_outs[s], row_vals[s], scale); break;
                case Bits::B32: write_quantized_vals<int32_t>(x_outs[s], row_vals[s], scale); break;
                case Bits::B99: write_quantized_vals<float>  (x_outs[s], row_vals[s], scale); break;
                default: throw std::runtime_error("Unsupported bits value");
            }
        }

        // Delta-encode row indices, matching the scheme in LDZipCompressor::write_i()
        using T = int16_t;
        static constexpr T DELTA_SENTINEL = std::numeric_limits<T>::max();

        std::vector<T> deltas(nnz);

        // First delta: global_col - row[0]  (stored as uint64 via int64 cast, same as write_i)
        uint64_t delta = static_cast<uint64_t>(
            static_cast<int64_t>(global_col) - static_cast<int64_t>(global_rows[0]));
        if (delta >= static_cast<uint64_t>(DELTA_SENTINEL)) {
            coo_out.push(0, global_col, static_cast<uint32_t>(delta));
            deltas[0] = DELTA_SENTINEL;
        } else {
            deltas[0] = static_cast<T>(delta);
        }

        // Subsequent deltas: row[k] - row[k-1]
        for (size_t idx = 1; idx < nnz; ++idx) {
            delta = static_cast<uint64_t>(
                static_cast<int64_t>(global_rows[idx]) -
                static_cast<int64_t>(global_rows[idx - 1]));
            if (delta >= static_cast<uint64_t>(DELTA_SENTINEL)) {
                coo_out.push(idx, global_col, static_cast<uint32_t>(delta));
                deltas[idx] = DELTA_SENTINEL;
            } else {
                deltas[idx] = static_cast<T>(delta);
            }
        }

        i_out.write(reinterpret_cast<const char *>(deltas.data()),
                    static_cast<std::streamsize>(nnz * sizeof(T)));
    }

    current_nnz += nnz;
    p_vec[global_col + 1] = current_nnz;
}

// ---------------------------------------------------------------------------
// concat_ldzip_overlapping
// ---------------------------------------------------------------------------
// Concatenates multiple .ldzip chunks that may have overlapping variant ranges.
// Overlap is auto-detected by comparing variant IDs in consecutive vars files.
// For each column in the overlapping region: entries with row <= column come
// from the earlier chunk (above-diagonal), entries with row > column come from
// the later chunk (below-diagonal).
// Only supports FULL matrix format (produced by plinkTabular compression).
void concat_ldzip_overlapping(const std::vector<std::string> &prefixes,
                              const std::string &out_prefix) {

    size_t n_chunks = prefixes.size();
    if (n_chunks == 0)
        throw std::runtime_error("No input chunks provided");

    // -----------------------------------------------------------------------
    // Step 1: Read vars files and detect pairwise overlaps
    // -----------------------------------------------------------------------
    std::vector<std::vector<std::string>> all_var_lines(n_chunks);
    for (size_t ci = 0; ci < n_chunks; ++ci)
        all_var_lines[ci] = read_var_lines(prefixes[ci] + ".vars.txt");

    std::vector<uint32_t> overlap_sizes(n_chunks > 1 ? n_chunks - 1 : 0, 0);
    for (size_t ci = 0; ci + 1 < n_chunks; ++ci) {
        overlap_sizes[ci] =
            compute_overlap_size(all_var_lines[ci], all_var_lines[ci + 1]);
        std::cout << " Overlap between chunk " << ci << " and chunk " << (ci + 1)
                  << ": " << overlap_sizes[ci] << " variants\n";
    }

    // -----------------------------------------------------------------------
    // Step 2: Compute global offsets (start index in the merged variant list
    //         for each chunk's first variant)
    // -----------------------------------------------------------------------
    std::vector<uint32_t> global_offsets(n_chunks, 0);
    for (size_t ci = 1; ci < n_chunks; ++ci) {
        global_offsets[ci] = global_offsets[ci - 1] +
            static_cast<uint32_t>(all_var_lines[ci - 1].size()) -
            overlap_sizes[ci - 1];
    }
    uint32_t total_n = global_offsets[n_chunks - 1] +
                       static_cast<uint32_t>(all_var_lines[n_chunks - 1].size());

    // -----------------------------------------------------------------------
    // Step 3: Open all input matrices and validate compatibility
    // -----------------------------------------------------------------------
    std::vector<LDZipMatrix> inputs;
    inputs.reserve(n_chunks);
    for (auto &p : prefixes)
        inputs.emplace_back(p);

    Bits bits = inputs[0].bitsEnum();
    MatrixFormat fmt = inputs[0].format();
    std::vector<Stat> stats = inputs[0].stats_available();

    if (fmt != MatrixFormat::FULL)
        throw std::runtime_error(
            "Overlapping concat is only supported for FULL matrix format "
            "(produced by plinkTabular compression).");

    for (size_t ci = 1; ci < n_chunks; ++ci) {
        if (inputs[ci].bitsEnum() != bits)
            throw std::runtime_error("Chunk " + std::to_string(ci) +
                                     " has a different bit width than chunk 0.");
        if (inputs[ci].format() != fmt)
            throw std::runtime_error("Chunk " + std::to_string(ci) +
                                     " has a different matrix format than chunk 0.");
        // Verify stats match
        auto other_stats = inputs[ci].stats_available();
        if (other_stats.size() != stats.size())
            throw std::runtime_error("Chunk " + std::to_string(ci) +
                                     " has different statistics than chunk 0.");
        for (size_t si = 0; si < stats.size(); ++si)
            if (other_stats[si] != stats[si])
                throw std::runtime_error("Chunk " + std::to_string(ci) +
                                         " has different statistics than chunk 0.");
    }

    // -----------------------------------------------------------------------
    // Step 4: Open output files
    // -----------------------------------------------------------------------
    std::string coo_path = out_prefix + ".io.bin";
    COO coo_out(coo_path.c_str(), 'w');

    std::ofstream i_out(out_prefix + ".i.bin",
                        std::ios::binary | std::ios::trunc);
    if (!i_out)
        throw std::runtime_error("Cannot create i file: " + out_prefix + ".i.bin");

    EnumArray<std::ofstream, Stat> x_outs;
    for (Stat s : stats) {
        std::string xp = out_prefix + ".x." + stat_to_string(s) + ".bin";
        x_outs[s].open(xp, std::ios::binary | std::ios::trunc);
        if (!x_outs[s])
            throw std::runtime_error("Cannot create x file: " + xp);
    }

    // -----------------------------------------------------------------------
    // Step 5: Process all global columns
    // -----------------------------------------------------------------------
    uint64_t current_nnz = 0;
    std::vector<uint64_t> p_vec(total_n + 1, 0);
    uint32_t global_col = 0;

    for (size_t ci = 0; ci < n_chunks; ++ci) {
        uint32_t n_cols   = static_cast<uint32_t>(all_var_lines[ci].size());
        uint32_t ov_prev  = (ci > 0)             ? overlap_sizes[ci - 1] : 0;
        uint32_t ov_next  = (ci + 1 < n_chunks)  ? overlap_sizes[ci]     : 0;

        // Non-overlap columns unique to this chunk (skip the leading columns
        // already emitted as part of the previous overlap pass, and the trailing
        // columns that will be emitted as part of the next overlap pass).
        uint32_t non_ov_start = ov_prev;
        uint32_t non_ov_end   = n_cols - ov_next;  // exclusive

        // --- Non-overlap columns ---
        for (uint32_t j = non_ov_start; j < non_ov_end; ++j) {
            auto local_rows = inputs[ci].get_i(j);
            EnumArray<std::vector<float>, Stat> row_vals;
            for (Stat s : stats)
                row_vals[s] = inputs[ci].get_x(j, s);

            // Translate local row indices to global
            std::vector<uint32_t> global_rows(local_rows.size());
            for (size_t k = 0; k < local_rows.size(); ++k)
                global_rows[k] = global_offsets[ci] + local_rows[k];

            write_column_to_files(global_col, global_rows, row_vals,
                                  current_nnz, p_vec,
                                  i_out, coo_out, x_outs, bits, stats);
            ++global_col;
        }

        // --- Overlap columns with the next chunk ---
        if (ci + 1 < n_chunks) {
            size_t ci_next = ci + 1;
            for (uint32_t ov = 0; ov < ov_next; ++ov) {
                // j1: local column index in chunk ci  (tail of chunk ci)
                // j2: local column index in chunk ci+1 (head of chunk ci+1)
                uint32_t j1 = non_ov_end + ov;   // = n_cols - ov_next + ov
                uint32_t j2 = ov;

                auto local_rows1 = inputs[ci].get_i(j1);
                auto local_rows2 = inputs[ci_next].get_i(j2);

                EnumArray<std::vector<float>, Stat> rv1, rv2;
                for (Stat s : stats) {
                    rv1[s] = inputs[ci].get_x(j1, s);
                    rv2[s] = inputs[ci_next].get_x(j2, s);
                }

                // Merge rule:
                //   From chunk ci  : rows with local_row <= j1
                //                    (global row <= global_col → above-diagonal)
                //   From chunk ci+1: rows with local_row >  j2
                //                    (global row >  global_col → below-diagonal)
                std::vector<uint32_t> merged_rows;
                EnumArray<std::vector<float>, Stat> merged_vals;

                for (size_t k = 0; k < local_rows1.size(); ++k) {
                    if (local_rows1[k] <= j1) {
                        merged_rows.push_back(global_offsets[ci] + local_rows1[k]);
                        for (Stat s : stats)
                            merged_vals[s].push_back(rv1[s][k]);
                    }
                }
                for (size_t k = 0; k < local_rows2.size(); ++k) {
                    if (local_rows2[k] > j2) {
                        merged_rows.push_back(global_offsets[ci_next] + local_rows2[k]);
                        for (Stat s : stats)
                            merged_vals[s].push_back(rv2[s][k]);
                    }
                }

                // merged_rows is already sorted ascending:
                // chunk ci rows are <= global_col; chunk ci+1 rows are > global_col

                write_column_to_files(global_col, merged_rows, merged_vals,
                                      current_nnz, p_vec,
                                      i_out, coo_out, x_outs, bits, stats);
                ++global_col;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Step 6: Write p-file
    // -----------------------------------------------------------------------
    {
        std::ofstream p_out(out_prefix + ".p.bin",
                            std::ios::binary | std::ios::trunc);
        if (!p_out)
            throw std::runtime_error("Cannot create p file: " + out_prefix + ".p.bin");
        p_out.write(reinterpret_cast<const char *>(p_vec.data()),
                    static_cast<std::streamsize>((total_n + 1) * sizeof(uint64_t)));
    }

    // -----------------------------------------------------------------------
    // Step 7: Write metadata JSON
    // -----------------------------------------------------------------------
    MetaInfo meta(total_n, total_n, current_nnz,
                  bits_to_int(bits), fmt, LDZipMatrix::DEFAULT_VERSION);
    for (Stat s : stats)
        meta.has_stat[s] = true;
    write_metadata_json(out_prefix + ".meta.json", meta);

    // -----------------------------------------------------------------------
    // Step 8: Merge vars file (deduplicated)
    // -----------------------------------------------------------------------
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