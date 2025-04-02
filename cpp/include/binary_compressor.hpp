#pragma once
#include <ldzipcompressor.hpp>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace ldzip {
    
    void compress_binary_matrix(const std::string& binary_path,
                                const std::string& vars_path,
                                const std::string& output_prefix,
                                int bits = 8,
                                float min = 1e-4,
                                const std::string& format_str = "full",
                                const Stat& stat = Stat::PHASED_R){
        if (bits != 8 && bits != 16 && bits != 32 && bits != 99) {
            throw std::runtime_error("Only bit sizes 8, 16, 32, or 99 are supported.");
        }
        
        MatrixFormat format = parse_format(format_str);
        std::ifstream in(binary_path, std::ios::binary);
        if (!in.is_open()) throw std::runtime_error("Failed to open input file: " + binary_path);

        // Step 1: Find first non-all-NA column and remove them
        size_t N = snp_util::count_lines(vars_path);
        std::vector<float> buffer(N);
        size_t i = 0;
        bool found = false;
        while (i < N && !found) {
            in.read(reinterpret_cast<char*>(&buffer[0]), N * sizeof(float));
            for (size_t j = 0; j < N; ++j) {
                if (!std::isnan(buffer[j])) {
                    found = true;
                    break;
                }
            }
            if (!found) ++i;
        }
        in.close();
        if (!found) throw std::runtime_error("All columns are empty/NA in "+binary_path);

        std::vector<size_t> keep;
        for (size_t j = 0; j < N; ++j) if (!std::isnan(buffer[j])) keep.push_back(j);
        size_t K = keep.size();


        std::cout << " found [" << (N) << "] variants in R2 matrix\n";
        std::cout << " found [" << (N - K) << "] bad variants and removed them \n";
        std::cout << " compressing to [" << K << " x " << K << "]\n";

        // Step 2: Read each column and compress them
        // LDZipMatrix ld(K, K, format, stat, parse_bits(bits), output_prefix, true);
        LDZipCompressor compressor(K, K, format, stat, parse_bits(bits), output_prefix, LDZipCompressor::Mode::ColumnStream);
        in.open(binary_path, std::ios::binary);

        for (size_t col_pos = 0; col_pos < K; ++col_pos) {
            if(col_pos%1000==0)
                std::cout << " compressing column [" << col_pos << "]\n";

            in.seekg(keep[col_pos] * N * sizeof(float), std::ios::beg);
            in.read(reinterpret_cast<char*>(&buffer[0]), N * sizeof(float));
            compressor.push_column(col_pos, buffer, keep, min, stat);
        }

        compressor.stream_close();
        in.close();
        snp_util::write_snp_subset(vars_path, output_prefix + ".vars.txt", keep);
    };

}
 
