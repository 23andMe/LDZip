#include "metadata.hpp"
#include "ldzipmatrix.hpp"
#include "ldzipcompressor.hpp"
#include "ldzipconcatenator.hpp"
#include "snp_util.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>
#include <stdexcept>

namespace ldzip {


void concat_ldzip(const std::vector<std::string> &prefixes,
                        const std::string &out_prefix,
                        bool overlapping) {

    std::vector<std::string> var_files;
    for (const auto& pref : prefixes)
        var_files.push_back(pref + ".vars.txt");

    // Build variant boundary information
    // If overlapping: parses files and validates overlaps
    // If non-overlapping: treats chunks as independent
    OverlapVariantInfo ov = read_overlapping_variant_order(var_files, overlapping);

    size_t total_rows = ov.total_variants;

    LDZipMatrix in(prefixes[0]);
    LDZipConcatenator concator(total_rows, total_rows, in.format(), in.stats_available(), in.bitsEnum(), out_prefix);

    // Process columns chunk by chunk
    for (size_t i = 0; i < prefixes.size(); i++) {
        std::cout << " Processing chunk " << i << " - "<< prefixes[i] << "\n";
        LDZipMatrix in(prefixes[i]);

        // --- Exclusive columns: raw binary copy (no decode/re-encode) ---
        concator.process_exclusive_columns(in, ov.chunks[i]);

        // --- Overlap columns with next chunk ---
        if (overlapping && i + 1 < prefixes.size()) {
            LDZipMatrix next_mat(prefixes[i+1]);
            concator.process_overlap_columns(in, next_mat, ov.chunks[i], ov.chunks[i + 1]);
        }
    }

    // Close streams and write metadata
    concator.close();

    // Write merged vars.txt
    // For overlapping: skip first_overlap_end variants (already written by previous chunk)
    // For non-overlapping: first_overlap_end=0, so writes all variants
    {
        std::string out_vars = out_prefix + ".vars.txt";
        std::ofstream out(out_vars);
        if (!out) throw std::runtime_error("Cannot open: " + out_vars);

        for (size_t i = 0; i < prefixes.size(); i++) {
            std::ifstream in(var_files[i]);
            if (!in) throw std::runtime_error("Cannot open: " + var_files[i]);
            const auto& chunk = ov.chunks[i];
            std::string line;
            size_t line_idx = 0;
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '#') {
                    if (i == 0) out << line << '\n';
                    continue;
                }
                // Write line if past the overlap with previous chunk
                if (line_idx >= chunk.first_overlap_end)
                    out << line << '\n';
                ++line_idx;
            }
        }
    }
}


}