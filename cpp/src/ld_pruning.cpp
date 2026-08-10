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
    uint32_t global_idx;

    static VariantInfo parse(const std::string& chrom_str, const std::string& pos_str, uint32_t idx) {
        VariantInfo v;
        v.global_idx = idx;
        v.pos = std::stoull(pos_str);

        std::string chrom = chrom_str;
        if (chrom.substr(0, 3) == "chr") {
            chrom = chrom.substr(3);
        }

        if (chrom == "X") {
            v.chrom = 23;
        } else if (chrom == "Y") {
            v.chrom = 24;
        } else if (chrom == "XY") {
            v.chrom = 25;
        } else if (chrom == "MT" || chrom == "M") {
            v.chrom = 26;
        } else {
            v.chrom = std::stoul(chrom);
        }

        return v;
    }
};

static size_t process_chromosome(
    std::vector<VariantInfo>& chrom_variants,
    uint32_t current_chrom,
    LDZipMatrix& matrix,
    double threshold,
    Stat stat,
    uint64_t window_bp,
    std::ofstream& prune_in,
    std::ofstream& prune_out
) {
    if (chrom_variants.empty()) return 0;

    std::vector<bool> removed(chrom_variants.size(), false);
    size_t removed_count = 0;

    for (uint32_t i = 0; i < chrom_variants.size(); ++i) {
        if (removed[i]) continue;

        uint32_t global_idx = chrom_variants[i].global_idx;
        const auto& neighbors = matrix.get_neighbors(global_idx, threshold, stat);

        for (uint32_t neighbor_idx : neighbors) {
            if (neighbor_idx <= global_idx) continue;

            for (uint32_t j = i + 1; j < chrom_variants.size(); ++j) {
                if (chrom_variants[j].global_idx != neighbor_idx) continue;
                if (removed[j]) break;

                uint64_t dist_bp = chrom_variants[j].pos - chrom_variants[i].pos;
                if (dist_bp > window_bp) break;

                removed[j] = true;
                removed_count++;
                break;
            }
        }
    }

    for (uint32_t i = 0; i < chrom_variants.size(); ++i) {
        if (removed[i]) {
            prune_out << chrom_variants[i].global_idx << "\n";
        } else {
            prune_in << chrom_variants[i].global_idx << "\n";
        }
    }

    std::cout << "Chromosome " << current_chrom << ": kept " << (chrom_variants.size() - removed_count) << ", removed " << removed_count << std::endl;

    chrom_variants.clear();
    return removed_count;
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

    bool use_full_preload = (matrix.get_total_file_size() <= 8ULL * 1024 * 1024 * 1024);
    if (use_full_preload) {
        matrix.preload_data();
    }

    const uint64_t window_bp = static_cast<uint64_t>(window_kb) * 1000;

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

    size_t total_variants = 0;
    size_t total_removed = 0;
    std::string line;
    std::vector<VariantInfo> chrom_variants;
    uint32_t current_chrom = 0;
    bool first_variant = true;
    uint32_t variant_idx = 0;

    while (std::getline(vars_in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string chrom_str, pos_str;
        if (!(ss >> chrom_str >> pos_str)) {
            throw std::runtime_error("Invalid variant file format");
        }

        VariantInfo v = VariantInfo::parse(chrom_str, pos_str, variant_idx);

        if (first_variant) {
            current_chrom = v.chrom;
            first_variant = false;
        } else if (v.chrom != current_chrom) {
            if (!chrom_variants.empty()) {
                std::cout << "Processing chromosome " << current_chrom << " (" << chrom_variants.size() << " variants)..." << std::endl;
                if (!use_full_preload) {
                    matrix.preload_segment(chrom_variants.front().global_idx, chrom_variants.back().global_idx);
                }
                total_removed += process_chromosome(chrom_variants, current_chrom, matrix, threshold, stat, window_bp, prune_in, prune_out);
            }
            current_chrom = v.chrom;
        }

        chrom_variants.push_back(v);
        total_variants++;
        variant_idx++;
    }

    if (!chrom_variants.empty()) {
        std::cout << "Processing chromosome " << current_chrom << " (" << chrom_variants.size() << " variants)..." << std::endl;
        if (!use_full_preload) {
            matrix.preload_segment(chrom_variants.front().global_idx, chrom_variants.back().global_idx);
        }
        total_removed += process_chromosome(chrom_variants, current_chrom, matrix, threshold, stat, window_bp, prune_in, prune_out);
    }

    vars_in.close();

    prune_in.close();
    prune_out.close();

    if (total_variants != matrix.nrows()) {
        throw std::runtime_error("Variant count mismatch");
    }

    std::cout << "Pruning complete: kept " << (total_variants - total_removed) << ", removed " << total_removed << std::endl;
    std::cout << "Output written to " << output_prefix << ".prune.in and .prune.out" << std::endl;
}

} // namespace ldzip
