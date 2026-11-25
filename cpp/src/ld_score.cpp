#include "ld_score.hpp"
#include "ldzipmatrix.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
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

void ld_score(
    const std::string& input_prefix,
    const std::string& output_file,
    size_t window_kb,
    float threshold,
    const std::string& stat_str
) {
    std::cout << "Loading matrix from " << input_prefix << "..." << std::endl;
    LDZipMatrix matrix(input_prefix);

    Stat stat = parse_stat(stat_str);
    if (!matrix.has_stat(stat)) {
        throw std::runtime_error("Statistic " + stat_str + " not available in matrix");
    }

    bool need_square = false;
    if (stat == Stat::PHASED_R || stat == Stat::UNPHASED_R) {
        need_square = true;
    } else if (stat == Stat::PHASED_R2 || stat == Stat::UNPHASED_R2) {
        need_square = false;
    } else {
        throw std::runtime_error("LD score only supports R or R2 statistics, got: " + stat_str);
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

    const uint64_t window_bp = static_cast<uint64_t>(window_kb) * 1000;

    std::ofstream out(output_file);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + output_file);
    }

    const size_t buffer_size = 10 * 1024 * 1024;
    std::vector<char> out_buffer(buffer_size);
    out.rdbuf()->pubsetbuf(out_buffer.data(), buffer_size);

    out << "variant\tld_score\n";

    for (uint32_t i = 0; i < variants.size(); ++i) {
        const auto& i_buf = matrix.get_i(i);
        const auto& x_buf = matrix.get_x(i, stat);

        double ld_score_val = 0.0;

        for (size_t k = 0; k < i_buf.size(); ++k) {
            uint32_t j = i_buf[k];
            float ld_val = x_buf[k];

            if (i == j) continue;
            if (std::abs(ld_val) < threshold) continue;
            if (variants[i].chrom != variants[j].chrom) continue;

            uint64_t dist_bp = (variants[j].pos >= variants[i].pos) ? (variants[j].pos - variants[i].pos) : (variants[i].pos - variants[j].pos);
            if (dist_bp > window_bp) continue;

            double r2 = need_square ? (ld_val * ld_val) : ld_val;
            ld_score_val += r2;
        }

        out << i << "\t" << ld_score_val << "\n";

        if ((i + 1) % 10000 == 0) {
            std::cout << "Processed " << (i + 1) << " / " << variants.size() << " variants" << std::endl;
        }
    }

    out.close();

    std::cout << "LD scores written to " << output_file << std::endl;
}

} // namespace ldzip
