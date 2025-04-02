#pragma once
#include "ldzipmatrix.hpp"
#include <string>
#include <fstream>
#include <stdexcept>

namespace ldzip {

    void decompress_ldzip(const std::string& prefix, const std::string& output_prefix, const std::string& type) {
        // Read metadata
        LDZipMatrix in(prefix);
        if(type=="tabular")
            in.readVariants(prefix + ".vars.txt");
        in.decompress(output_prefix, type);
        
        
        // Copy the SNP file
        std::ifstream snp_in(prefix + ".vars.txt");
        std::string snp_out_path = (type == "tabular") ? output_prefix + ".vars.txt" : output_prefix + ".bin.vars";
        std::ofstream snp_out(snp_out_path);
        if (!snp_in || !snp_out)
            throw std::runtime_error("Failed to open SNP file for copy");
        snp_out << snp_in.rdbuf();
    };
    
}
