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

static uint32_t parse_chrom(const std::string& chrom_str) {
    std::string chrom = chrom_str;
    if (chrom.substr(0, 3) == "chr") {
        chrom = chrom.substr(3);
    }

    if (chrom == "X") return 23;
    if (chrom == "Y") return 24;
    if (chrom == "XY") return 25;
    if (chrom == "MT" || chrom == "M") return 26;
    return std::stoul(chrom);
}

void ld_pruning(
    const std::string& input_prefix,
    const std::string& output_prefix,
    double threshold,
    size_t window_kb,
    const std::string& stat_str
) {
    std::cout << "Loading matrix from " << input_prefix << "..." << std::endl;
    LDZipMatrix matrix(input_prefix);

    Stat stat = parse_stat(stat_str);
    if (!matrix.has_stat(stat)) {
        throw std::runtime_error("Statistic " + stat_str + " not available in matrix");
    }

    std::string vars_file = input_prefix + ".vars.txt";
    std::ifstream vars_in(vars_file);
    if (!vars_in) {
        throw std::runtime_error("Cannot open variant file: " + vars_file);
    }

    std::vector<VariantInfo> variants;
    variants.reserve(matrix.nrows());

    std::cout << "Loading variants..." << std::endl;
    std::string line;
    while (std::getline(vars_in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string chrom_str, pos_str;
        if (!(ss >> chrom_str >> pos_str)) {
            throw std::runtime_error("Invalid variant file format");
        }

        VariantInfo v;
        v.chrom = parse_chrom(chrom_str);
        v.pos = std::stoull(pos_str);
        variants.push_back(v);
    }
    vars_in.close();

    if (variants.size() != matrix.nrows()) {
        throw std::runtime_error("Variant count mismatch");
    }

    std::cout << "Variants: " << variants.size() << std::endl;

    if (matrix.get_total_file_size() <= 8ULL * 1024 * 1024 * 1024) {
        matrix.preload_data();
    }

    std::vector<bool> removed(variants.size(), false);
    size_t removed_count = 0;

    const uint64_t window_bp = static_cast<uint64_t>(window_kb) * 1000;

    for (uint32_t i = 0; i < variants.size(); ++i) {
        if (removed[i]) continue;

        const auto& neighbors = matrix.get_neighbors(i, threshold, stat);

        for (uint32_t j : neighbors) {
            if (j <= i) continue;
            if (removed[j]) continue;
            if (variants[i].chrom != variants[j].chrom) {
                throw std::runtime_error("Cross-chromosome LD should not exist");
            }

            uint64_t dist_bp = (variants[j].pos >= variants[i].pos) ? (variants[j].pos - variants[i].pos) : (variants[i].pos - variants[j].pos);
            if (dist_bp > window_bp) continue;

            removed[j] = true;
            removed_count++;
        }

        if ((i + 1) % 10000 == 0) {
            std::cout << "Processed " << (i + 1) << " / " << variants.size() << " variants, removed " << removed_count << std::endl;
        }
    }

    std::cout << "Pruning complete: kept " << (variants.size() - removed_count) << ", removed " << removed_count << std::endl;

    std::ofstream prune_in(output_prefix + ".prune.in");
    std::ofstream prune_out(output_prefix + ".prune.out");

    if (!prune_in || !prune_out) {
        throw std::runtime_error("Failed to open output files");
    }

    const size_t buffer_size = 10 * 1024 * 1024;
    std::vector<char> in_buffer(buffer_size);
    std::vector<char> out_buffer(buffer_size);
    prune_in.rdbuf()->pubsetbuf(in_buffer.data(), buffer_size);
    prune_out.rdbuf()->pubsetbuf(out_buffer.data(), buffer_size);

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
