#include "find_tag_variants.hpp"
#include "ldzipmatrix.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace ldzip {

void find_tag_variants(
    const std::string& input_prefix,
    const std::string& variant_file,
    const std::string& output_file,
    double threshold,
    const std::string& stat_str
) {
    // Load the matrix
    std::cout << "Loading matrix from " << input_prefix << "..." << std::endl;
    LDZipMatrix matrix(input_prefix);

    Stat stat = parse_stat(stat_str);

    // Decide preload strategy based on file size
    const size_t max_full_preload = 8ULL * 1024 * 1024 * 1024; // 8GB
    size_t total_size = matrix.get_total_file_size();
    bool use_full_preload = (total_size <= max_full_preload);

    if (use_full_preload) {
        matrix.preload_data();
    }

    // Open input and output files
    std::ifstream in(variant_file);
    if (!in) {
        throw std::runtime_error("Failed to open variant file: " + variant_file);
    }

    std::ofstream out(output_file);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + output_file);
    }

    // Set large buffer for output stream
    const size_t buffer_size = 10 * 1024 * 1024;
    std::vector<char> out_buffer(buffer_size);
    out.rdbuf()->pubsetbuf(out_buffer.data(), buffer_size);

    // Write header
    out << "variant\ttag_variant\n";

    // Process variants line by line
    std::string line;
    int line_no = 0;
    size_t total_tags = 0;
    size_t variants_processed = 0;
    uint32_t current_segment_end = 0;

    std::cout << "Processing variants from " << variant_file << "..." << std::endl;

    while (std::getline(in, line)) {
        ++line_no;

        // Parse index from line
        std::istringstream iss(line);
        int idx;
        char extra;

        if (!(iss >> idx) || (iss >> extra)) {
            throw std::runtime_error("Invalid line " + std::to_string(line_no) +
                                   " in variant file (must contain exactly one integer)");
        }
        if (idx < 0) {
            throw std::runtime_error("Negative index at line " + std::to_string(line_no));
        }
        if (static_cast<uint32_t>(idx) >= matrix.nrows()) {
            throw std::runtime_error("Index " + std::to_string(idx) +
                                   " at line " + std::to_string(line_no) +
                                   " exceeds matrix size (" + std::to_string(matrix.nrows()) + ")");
        }

        uint32_t var_idx = static_cast<uint32_t>(idx);
        ++variants_processed;

        // Segment-based preloading for large files
        if (!use_full_preload && var_idx > current_segment_end) {
            const size_t segment_size = 4ULL * 1024 * 1024 * 1024; // 4GB segments
            current_segment_end = matrix.get_segment_end(var_idx, segment_size);
            matrix.preload_segment(var_idx, current_segment_end);
        }

        // Progress reporting
        if (variants_processed % 1000 == 0) {
            std::cout << "Processed " << variants_processed << " variants, found "
                      << total_tags << " tags so far..." << std::endl;
        }

        // Get neighbors (tags) for this variant
        std::vector<uint32_t> tags = matrix.get_neighbors(var_idx, threshold, stat);

        // Write tags (excluding self)
        for (uint32_t tag_idx : tags) {
            if (tag_idx != var_idx) {
                out << var_idx << "\t" << tag_idx << "\n";
                ++total_tags;
            }
        }
    }

    in.close();
    out.close();

    std::cout << "Done! Found " << total_tags << " tag variants for "
              << variants_processed << " input variants" << std::endl;
    std::cout << "Output written to " << output_file << std::endl;
}

} // namespace ldzip
