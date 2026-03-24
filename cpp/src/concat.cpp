#include "metadata.hpp"
#include "ldzipmatrix.hpp"
#include "ldzipcompressor.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>
#include <stdexcept>

namespace ldzip {

template <typename T>
void copy_binary_file_with_offset(const std::string& in_file,
                                  const std::string& out_file,
                                  std::streamoff read_pos,
                                  T add_value,
                                  size_t byte_size,
                                  size_t n,
                                  size_t chunk_size = 1000000) {
    std::ifstream in(in_file, std::ios::binary);
    if (!in) throw std::runtime_error("Failed to open input file: " + in_file);
    std::cout<<" copying "<<in_file<<" to "<<out_file<<std::endl;
    if (read_pos > 0) {
        in.seekg(read_pos, std::ios::beg);
    }

    std::ofstream out(out_file, std::ios::binary | std::ios::app);
    if (!out) throw std::runtime_error("Failed to open output file: " + out_file);

    size_t remaining = n;
    if (n == 0) {
        std::ifstream size_in(in_file, std::ios::binary | std::ios::ate);
        std::streamoff bytes = size_in.tellg() - read_pos;
        if (bytes < 0 || bytes % byte_size != 0)
            throw std::runtime_error("Invalid file size for element type");

        remaining = bytes / byte_size;
    }

    while (remaining > 0) {
        size_t to_read = (chunk_size < remaining) ? chunk_size : remaining;
        std::vector<T> buffer(to_read);

        in.read(reinterpret_cast<char*>(buffer.data()), to_read * byte_size);

        if constexpr (std::is_same_v<T, ldzip::COO::Entry>) {
            for (auto& e : buffer) {
                e.j += add_value.j;
            }
        } else {
            if (add_value != T{}) {
                for (auto& v : buffer) {
                    v += add_value;
                }
            }
        }

        out.write(reinterpret_cast<char*>(buffer.data()), buffer.size() * byte_size);
        remaining -= to_read;
    }
}


void merge_vars_files(const std::vector<std::string> &prefixes, const std::string &out_file) {
    std::remove(out_file.c_str());
    std::ofstream out(out_file);
    if (!out) throw std::runtime_error("Cannot open output file: " + out_file);

    bool firstFile = true;

    for (const auto &prefix : prefixes) {
        std::ifstream in(prefix + ".vars.txt");
        if (!in) throw std::runtime_error("Cannot open input file: " + prefix + ".vars.txt");

        std::string line;

        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '#')
            {
                if (firstFile) out << line << '\n';
                continue;
            }
            out << line << '\n';
        }
        firstFile = false;
    }
}


void concat_ldzip(const std::vector<std::string> &prefixes,
                        const std::string &out_prefix) {

    uint32_t row_offset = 0, total_rows = 0;
    uint64_t nnz_offset = 0, total_nnz = 0;
    Bits bits;

    for (const auto &pref : prefixes) {
        LDZipMatrix in(pref);
        total_rows += static_cast<uint32_t>(in.nrows());
        total_nnz += static_cast<uint64_t>(in.nnz());
    }

    LDZipMatrix in(prefixes[0]);
    LDZipCompressor concator(total_rows, total_rows, in.format(), in.stats_available(), in.bitsEnum(), out_prefix, LDZipCompressor::Mode::ColumnStream);
    concator.m().set_nnz(total_nnz);

    // Initialize p_file with 0
    {
        std::ofstream con_p(concator.m().pFile(), std::ios::binary | std::ios::app);
        uint64_t zero = 0;
        con_p.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
    }

    for (const auto &pref : prefixes) {
        std::cout << " reading : " << pref << "\n";
        LDZipMatrix in(pref);
        size_t nrows = in.nrows();
        bits = in.bitsEnum();
        uint64_t nnz = in.nnz();

        copy_binary_file_with_offset<uint64_t>(in.pFile(), concator.m().pFile(), sizeof(uint64_t), nnz_offset, sizeof(uint64_t), nrows);
        copy_binary_file_with_offset<int16_t>(in.iFile(), concator.m().iFile(), 0, 0, sizeof(int16_t), nnz);
        copy_binary_file_with_offset<ldzip::COO::Entry>(in.IFile(), concator.m().IFile(), 0, ldzip::COO::Entry{0, row_offset, 0}, sizeof(ldzip::COO::Entry), 0);
        for (Stat s : All_Stats()) 
            if (in.has_stat(s)) {
                if(bits!=Bits::B99)
                    copy_binary_file_with_offset<int32_t>(in.xFile(s), concator.m().xFile(s), 0L, 0, bits_to_int(bits) / 8, nnz);
                else
                    copy_binary_file_with_offset<float>(in.xFile(s), concator.m().xFile(s), 0L, 0.0f, sizeof(float), nnz);
            }
        row_offset += nrows;
        nnz_offset += nnz;
    }

    write_metadata_json(concator.m().metaFile(), concator.m().metaInfo());
    merge_vars_files(prefixes, out_prefix + ".vars.txt");
}


}