#pragma once
#include "ldzipcompressor.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <unordered_map>


namespace ldzip {
    
void compress_tabular_file(const std::string& in_file, 
                                const std::string& vars_path,
                                const std::string& output_prefix,
                                int bits,
                                float min,
                                const std::string& format_str,
                                Stat& min_stat){
        if (bits != 8 && bits != 16 && bits != 32 && bits != 99) {
            throw std::runtime_error("Only bit sizes 8, 16, 32, or 99 are supported.");
        }

        MatrixFormat format = parse_format(format_str);
        std::ifstream in(in_file, std::ios::in);
        if (!in.is_open()) throw std::runtime_error("Failed to open input file: " + in_file);

        // Get hash of SNPs to get their row and col indices efficiently
        size_t N = snp_util::count_lines(vars_path);
        int noFields = snp_util::countFields(in_file);
        std::cout<<" Generating hash of SNPs"<<std::endl;
        std::unordered_map<std::string, uint32_t> snp_hash;
        snp_util::get_snp_hash(vars_path, snp_hash);

        // Initialize the compressor
        std::vector<std::string> fields = snp_util::getFields(in_file, true);
        auto [stats_available, col_to_stat] = snp_util::get_stats_from_fields(fields);
        LDZipCompressor compressor(N, N, format, stats_available, parse_bits(bits), output_prefix, LDZipCompressor::Mode::ValueStream);


        // Temporary variables to start reading data
        EnumArray<float, Stat> values(-999.0f); 
        uint32_t c, r; 
        std::string line;
        uint32_t start=0;
        uint64_t count=0;

        getline(in, line); // skip header

        while (getline(in, line)) {
            if (count++ % 1000000 == 0)
                std::cout << "Processed " << count-1 << " entries" << std::endl;
            std::stringstream ss(line);
            for (int i = 0; i < noFields; ++i) {
                getline(ss, fields[i], '\t');
            }

            if (snp_hash.count(fields[0]) == 0 || snp_hash.count(fields[3]) == 0)
                throw std::runtime_error("SNP not found in hash: " + fields[0] + " or " + fields[3]);
                
            // read in values
            c = snp_hash[fields[0]] - start;
            r = snp_hash[fields[3]] - start;
            for (auto& [col, stat] : col_to_stat) {
                values[stat] = std::stof(fields[col]);
            }

            // check values
            if (values[min_stat] == -999.0f) {
                throw std::runtime_error("Column " + stat_to_string(min_stat) + " not available, cannot use [--min_col "+ stat_to_string(min_stat)+"]");
            }
            if (std::abs(values[min_stat]) < min) continue;

            // push values
            compressor.push_value(r, c, values);
        }

        compressor.stream_close();
        in.close();
        snp_util::write_snp_subset(vars_path, output_prefix + ".vars.txt");

        
    };

}
