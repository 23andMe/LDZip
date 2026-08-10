#include "ld_pruning.hpp"
#include "ldzipmatrix.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

namespace ldzip {

struct VariantInfo {
    uint32_t chrom;
    uint64_t pos;
};

// LD pruning implementation
//
// This is a forward greedy algorithm: variants are processed in increasing index order.
// For each surviving variant i, all later variants j > i within the genomic window
// that have |LD(i,j)| >= threshold are removed.
//
// The first variant in each high-LD pair is retained, the second is removed.
// This is intentionally NOT identical to PLINK --indep-pairwise, which uses:
//   - Backward scanning within windows
//   - MAF-based tie-breaking (prefers keeping higher MAF variants)
//
// This implementation provides deterministic LD pruning but may produce
// different results than PLINK for the same threshold.
void ld_pruning(
    const std::string& input_prefix,
    const std::string& output_prefix,
    double threshold,
    size_t window_kb,
    const std::string& stat_str
) {
    // Load the matrix
    std::cout << "Loading matrix from " << input_prefix << "..." << std::endl;
    LDZipMatrix matrix(input_prefix);

    Stat stat = parse_stat(stat_str);
    if (!matrix.has_stat(stat)) {
        throw std::runtime_error("Statistic " + stat_str + " not available in matrix");
    }

    // Read variant positions from .vars.txt
    std::string vars_file = input_prefix + ".vars.txt";
    std::ifstream vars_in(vars_file);
    if (!vars_in) {
        throw std::runtime_error("Cannot open variant file: " + vars_file);
    }

    std::vector<VariantInfo> variants;
    variants.reserve(matrix.nrows());

    std::string line;
    while (std::getline(vars_in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string chrom_str, pos_str;
        if (!(ss >> chrom_str >> pos_str)) {
            throw std::runtime_error("Invalid variant file format");
        }

        VariantInfo v;
        v.chrom = std::stoul(chrom_str);
        v.pos = std::stoull(pos_str);
        variants.push_back(v);
    }
    vars_in.close();

    if (variants.size() != matrix.nrows()) {
        throw std::runtime_error("Variant count mismatch");
    }

    std::cout << "Variants: " << variants.size() << std::endl;
    std::cout << "Window: " << window_kb << " kb" << std::endl;
    std::cout << "Threshold: " << threshold << std::endl;
    std::cout << "Statistic: " << stat_str << std::endl;

    // Preload data if small enough
    const size_t max_full_preload = 8ULL * 1024 * 1024 * 1024; // 8GB
    size_t total_size = matrix.get_total_file_size();
    bool use_full_preload = (total_size <= max_full_preload);

    if (use_full_preload) {
        matrix.preload_data();
    }

    // Track removed variants
    std::vector<bool> removed(variants.size(), false);
    size_t removed_count = 0;

    // Compute window size in bp once
    const uint64_t window_bp = static_cast<uint64_t>(window_kb) * 1000;

    // Process each variant
    for (uint32_t i = 0; i < variants.size(); ++i) {
        if (removed[i]) continue;

        // Get neighbors within window that exceed threshold
        const auto& neighbors = matrix.get_neighbors(i, threshold, stat);

        for (uint32_t j : neighbors) {
            if (j <= i) continue;
            if (removed[j]) continue;
            if (variants[i].chrom != variants[j].chrom) continue;

            uint64_t dist_bp = (variants[j].pos >= variants[i].pos)
                ? (variants[j].pos - variants[i].pos)
                : (variants[i].pos - variants[j].pos);

            if (dist_bp > window_bp) continue;

            removed[j] = true;
            removed_count++;
        }

        // Progress reporting
        if ((i + 1) % 100 == 0) {
            std::cout << "Processed " << (i + 1) << " / " << variants.size()
                      << " variants, removed " << removed_count << std::endl;
        }
    }

    std::cout << "Pruning complete: kept " << (variants.size() - removed_count)
              << ", removed " << removed_count << std::endl;

    // Write output files
    std::ofstream prune_in(output_prefix + ".prune.in");
    std::ofstream prune_out(output_prefix + ".prune.out");

    if (!prune_in || !prune_out) {
        throw std::runtime_error("Failed to open output files");
    }

    for (size_t i = 0; i < variants.size(); ++i) {
        if (removed[i]) {
            prune_out << i << "\n";
        } else {
            prune_in << i << "\n";
        }
    }

    prune_in.close();
    prune_out.close();

    std::cout << "Output written to " << output_prefix << ".prune.in and .prune.out" << std::endl;
}

} // namespace ldzip
