#pragma once

#include "metadata.hpp"
#include "snp_util.hpp"
#include "coo.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
static_assert(sizeof(float) == 4,   "float must be 4 bytes");
static_assert(sizeof(double) == 8,  "double must be 8 bytes");

namespace ldzip {

    class LDZipMatrix {
        public:
            static constexpr const char* DEFAULT_VERSION = "2.1";

            // --- Constructors ---
            LDZipMatrix     ();
            LDZipMatrix     (size_t nrows, size_t ncols, MatrixFormat format, Stat stat, Bits bits, const std::string& prefix);
            LDZipMatrix     (size_t nrows, size_t ncols, MatrixFormat format, const std::vector<Stat>& stats, Bits bits, const std::string& prefix);
            LDZipMatrix     (const std::string& prefix);

            // --- Metadata/Size Access ---
            std::string     version()       const noexcept { return version_; }
            size_t          nrows()         const noexcept { return nrows_; }
            size_t          ncols()         const noexcept { return ncols_; }
            uint64_t        nnz()           const noexcept {return nnz_; }
            void            set_nnz         (uint64_t nnz){nnz_=nnz; }
            bool            empty()         const noexcept { return nrows_ == 0 || ncols_ == 0; }

            // --- Stats and Metadata
            bool            has_stat(Stat s)    const noexcept { return has_stat_[s]; }
            bool            has_stat(const std::string& stat)    
                                                const noexcept { return has_stat(parse_stat(stat)); }
            std::vector<Stat> stats_available() const noexcept { return stats_available_; }

            Bits            bitsEnum()          const noexcept { return bits_; }
            int             bits()              const noexcept { return static_cast<int>(bits_); }
            MatrixFormat    format()            const noexcept { return format_; }
            std::string     formatStr()         const noexcept { return format_to_string(format_); }
            MetaInfo        metaInfo()          const;
            
            // --- Variant Metadata
            const std::vector<Variant>&     variants()      const { return variants_; }
            void                            readVariants    (const std::string& snp_file);

            // -- File Access ---
            bool            checkFiles()        const;
            bool            checkStatFiles()    const;
            bool            checkOverflowFiles()const;
            std::string     prefix()            const noexcept { return file_prefix_; }
            std::string     metaFile()          const noexcept { return file_prefix_ + fileSuffix(FileType::METADATA); }
            std::string     pFile()             const noexcept { return file_prefix_ + fileSuffix(FileType::P_VECTOR); }
            std::string     iFile()             const noexcept { return file_prefix_ + fileSuffix(FileType::I_VECTOR); }
            std::string     IFile()             const noexcept { return file_prefix_ + fileSuffix(FileType::I_OVERFLOW_VECTOR); }
            std::string     xFile(Stat s)       const noexcept { return file_prefix_ + xSuffix(s); }

            // --- Low Level vector access ---
            const std::vector<uint64_t>&        get_p()                             const;
            uint64_t                            get_p(size_t idx)                   const;
            std::vector<uint32_t>               get_i(uint32_t column)              const;
            std::vector<float>                  get_x(uint32_t column, Stat stat)   const;

            // --- Element Access ---
            float getValue(size_t row, size_t col, const Stat& stat) const;
            float getValue(size_t row, size_t col, const std::string& stat) const { return getValue(row, col, parse_stat(stat)); };

            // --- Column/Row Access ---
            void getColumn(uint32_t column, std::vector<float>& value, const Stat& stat) const;
            void getColumn(size_t col, size_t row_start, size_t row_end, double* values, const Stat& stat) const;
            void getColumn(size_t col, const int* rows, size_t nrows_sel, double* values, const Stat& stat) const;
            void getColumn(uint32_t column, std::vector<float>& value, const std::string& stat) const {
                getColumn(column, value, parse_stat(stat));
            };
            void getColumn(size_t col, size_t row_start, size_t row_end, double* values, const std::string& stat) const {
                getColumn(col, row_start, row_end, values, parse_stat(stat));
            };
            void getColumn(size_t col, const int* rows, size_t nrows_sel, double* values, const std::string& stat) const {
                getColumn(col, rows, nrows_sel, values, parse_stat(stat));
            };
            void getColumnRaw(size_t column, std::vector<size_t> rows, std::vector<size_t> &new_i, std::vector<float> &new_x, const Stat &stat) const;

            // --- Sub Matrix Access
            void getSubMatrix   (size_t col_start,      size_t col_end,
                                size_t row_start,       size_t row_end,
                                double* values,         const Stat& stat) const;
            void getSubMatrix(  const int* cols,        size_t ncols_sel,
                                const int* rows,        size_t nrows_sel,
                                double* values,         const Stat& stat) const;
            void getSubMatrix(  size_t col_start,       size_t col_end,
                                size_t row_start,       size_t row_end,
                                double* values,         const std::string& stat) const {
                getSubMatrix(col_start, col_end, row_start, row_end, values, parse_stat(stat));
            };
            void getSubMatrix(  const int* cols,        size_t ncols_sel,
                                const int* rows,        size_t nrows_sel,
                                double* values,         const std::string& stat) const{
                getSubMatrix(cols, ncols_sel, rows, nrows_sel, values, parse_stat(stat));
            };
            

            // --- Pairwise Access

            void getPairwise(const int* cols, const int* rows, size_t npairs, double* values, const Stat& stat) const;
            void getPairwise(const int* cols, const int* rows, size_t npairs, double* values, const std::string& stat) const {
                getPairwise(cols, rows, npairs, values, parse_stat(stat));
            };

            // --- Get Neighbors ---
            std::vector<uint32_t> get_neighbors(uint32_t column, double abs_threshold, const Stat& stat) const;
            std::vector<uint32_t> get_neighbors(uint32_t column, double abs_threshold, const std::string& stat) const {
                return get_neighbors(column, abs_threshold, parse_stat(stat));
            };

            // --- Compress and Decompress ---
            void    decompress              (const std::string& out_prefix, const std::string& type);
            void    writeColumnBinary       (std::ofstream& out, const std::vector<float>& value);
            void    writeColumnTabular      (std::ofstream& out, size_t column);

            // --- COO stream management for concatenation ---
            void    close_I()               { I_.close(); }
            void    reopen_I_append()       { I_.open_append(IFile().c_str()); }

        private:

            enum class FileType { I_VECTOR, I_OVERFLOW_VECTOR, P_VECTOR, METADATA };
            static constexpr const char* fileSuffix(FileType type) {
                switch (type) {
                    case FileType::I_VECTOR:            return ".i.bin";
                    case FileType::I_OVERFLOW_VECTOR:   return ".io.bin";
                    case FileType::P_VECTOR:            return ".p.bin";
                    case FileType::METADATA:            return ".meta.json";
                }
                return "";
            }

            static std::string xSuffix(Stat s) {
                return ".x." + stat_to_string(s) + ".bin";
            }

            // --- Private dimensions and metadata
            std::string version_;
            size_t nrows_{0};
            size_t ncols_{0};
            uint64_t nnz_{0};
            Bits bits_{Bits::B8};
            MatrixFormat format_{MatrixFormat::UPPER};
            mutable EnumArray<bool, Stat> has_stat_{};
            std::vector<Stat> stats_available_{};
            std::vector<Variant> variants_;

            //  --- Private file handlers 
            std::string file_prefix_{};
            mutable std::fstream p_stream_;
            mutable std::fstream i_stream_;
            mutable EnumArray<std::fstream, Stat> x_streams_;

            // --- Private Core vectors
            COO I_;             // COO style sparse matrix for i_ where Delta(i_) is outside uint16 range 
            mutable std::vector<uint64_t> p_;
            mutable std::vector<std::vector<uint32_t>> i_;
            mutable EnumArray<std::vector<std::vector<float>>, Stat> xs_;


            // --- Friend Class for Compressors            
            friend class LDZipCompressor;
            friend class LDZipConcatenator;


    };

} // namespace ldzip