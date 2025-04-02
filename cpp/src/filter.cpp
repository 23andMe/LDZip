#include "metadata.hpp"
#include "ldzipmatrix.hpp"
#include "ldzipcompressor.hpp"
#include <vector>
#include <string>

namespace ldzip {

void filter_ldzip(const std::string &in_prefix,
                        const std::string &out_prefix, 
                        std::vector<size_t> indices)
{
    std::cout << " filtering file : " << in_prefix << "\n";
    LDZipMatrix in(in_prefix);

    if(in.stats_available().size()>1)
        throw std::runtime_error("Only a single stat type can be filtered");
    if (indices.empty())
        throw std::runtime_error("Index list is empty");
    if (!std::is_sorted(indices.begin(), indices.end()))
        throw std::runtime_error("Indices must be sorted");
    if (indices.back() > in.nrows())
        throw std::runtime_error("Index exceeds maximum allowed value");

    Stat thisStat = in.stats_available()[0];

    LDZipCompressor compressor(indices.size(), 
                                indices.size(), 
                                in.format(), 
                                in.stats_available(), 
                                in.bitsEnum(), 
                                out_prefix, 
                                LDZipCompressor::Mode::ColumnStream);

    std::vector<size_t> new_i(indices.size());
    std::vector<float> new_x(indices.size());

    for (size_t idx = 0; idx < indices.size(); ++idx) {

        size_t col = indices[idx];
        in.getColumnRaw(col, indices, new_i, new_x, thisStat);
        compressor.push_column_raw(idx, new_x, new_i, thisStat);
    }

    compressor.stream_close();
    snp_util::write_snp_subset(in_prefix + ".vars.txt", out_prefix + ".vars.txt", indices);
}



}

