#include "ldzipmatrix.hpp"
#include <iostream>
#include <filesystem>
#include <cmath>
#include <sstream>
namespace ldzip {

// --- Helper Functions ---

// Check if an array is sorted in ascending order (allows duplicates)
template<typename T>
static void validate_sorted(const T* arr, size_t n, const char* name) {
    for (size_t i = 1; i < n; ++i) { if (arr[i] < arr[i-1]) { throw std::invalid_argument(std::string(name) + " array must be sorted in ascending order");}}
}

// Overload for std::vector
template<typename T>
static void validate_sorted(const std::vector<T>& vec, const char* name) {
    validate_sorted(vec.data(), vec.size(), name);
}

// --- Constructors ---
LDZipMatrix::LDZipMatrix() = default;

LDZipMatrix::LDZipMatrix(size_t nrows,
                         size_t ncols,
                         MatrixFormat format,
                         const std::vector<Stat>& stats,
                         Bits bits,
                         const std::string& prefix)
    : version_(DEFAULT_VERSION),
      nrows_(nrows),
      ncols_(ncols),
      nnz_(0),
      bits_(bits),
      format_(format),
      file_prefix_(prefix),
      I_(IFile().c_str(), 'w'){

    for (Stat s : stats) {
        has_stat_[s] = true;
        stats_available_.push_back(s);
    }

    p_.resize(ncols_ + 1);
    i_.resize(ncols_ + 1);
    for (Stat s : stats_available_) xs_[s].resize(ncols_ + 1);

}

LDZipMatrix::LDZipMatrix(size_t nrows,
                         size_t ncols,
                         MatrixFormat format,
                         Stat stat,
                         Bits bits,
                         const std::string& prefix)
    : LDZipMatrix(nrows, ncols, format, std::vector<Stat>{stat}, bits, prefix) {    
}


LDZipMatrix::LDZipMatrix(const std::string& prefix) :   file_prefix_(prefix),  
                                                        I_(IFile().c_str(), 'r') {
    checkFiles();

    MetaInfo meta = read_metadata_json(file_prefix_ + fileSuffix(FileType::METADATA));
    nrows_  = meta.rows;
    ncols_  = meta.cols;
    format_ = meta.format;
    bits_ = meta.bits;
    nnz_ = meta.nnz;
    has_stat_ = meta.has_stat;
    version_ = meta.version;
    checkOverflowFiles();
    I_.load();
    checkStatFiles();
    p_stream_.open(pFile(), std::ios::in | std::ios::binary);
    i_stream_.open(iFile(), std::ios::in | std::ios::binary);
    for (Stat s : All_Stats()) if (has_stat_[s]) {stats_available_.push_back(s); x_streams_[s].open(xFile(s), std::ios::in | std::ios::binary);}
}


// -- File Access ---
bool LDZipMatrix::checkFiles() const {
    for (auto type : {FileType::I_VECTOR, FileType::P_VECTOR, FileType::METADATA}) {
        std::string path = file_prefix_ + fileSuffix(type);
        if (!std::filesystem::exists(path)) {
            throw std::runtime_error("LDZipMatrix:: Missing file: " + path);
        }
    }
    return true;
}

bool LDZipMatrix::checkOverflowFiles() const {
    if (version_ != "1.1") {
        std::string path = file_prefix_ +
                           fileSuffix(FileType::I_OVERFLOW_VECTOR);
        if (!std::filesystem::exists(path)) {
            throw std::runtime_error("LDZipMatrix:: Missing file: " + path);
        }
    }
    return true;
}

bool LDZipMatrix::checkStatFiles() const {
    for (Stat s : stats_available_) 
        if (!std::filesystem::exists(xFile(s))) 
            throw std::runtime_error("LDZipMatrix:: Missing file: " + xFile(s));
    return true;
}

// --- Low Level vector access ---
const std::vector<uint64_t>& LDZipMatrix::get_p() const {
    if (!p_.empty()) return p_;

    p_.resize(ncols_ + 1);

    p_stream_.clear(); // clear EOF/fail flags if set
    p_stream_.seekg(0, std::ios::beg);
    p_stream_.read(reinterpret_cast<char*>(p_.data()), (nrows_ + 1) * sizeof(uint64_t));
    if (!p_stream_) throw std::runtime_error("Error reading p file: " + pFile());
    return p_;
} 

uint64_t LDZipMatrix::get_p(size_t idx) const {
    if (!p_.empty()) return p_[idx];
    uint64_t p_val;
    p_stream_.clear();
    p_stream_.seekg(idx * sizeof(uint64_t), std::ios::beg);
    p_stream_.read(reinterpret_cast<char*>(&p_val), sizeof(uint64_t));
    if (!p_stream_) throw std::runtime_error("Error reading p file: " + pFile());
    return p_val;
} 

std::vector<uint32_t> LDZipMatrix::get_i(uint32_t column) const {
    uint64_t start = get_p(column);
    uint64_t end = get_p(column+1);
    uint64_t nnz_column = end - start;
    std::vector<uint32_t> i_buf(nnz_column);
    if (nnz_column == 0) return i_buf;

    i_stream_.clear();

    if (version_ == "1.1") {
        i_stream_.seekg(start * sizeof(uint32_t), std::ios::beg);
        i_stream_.read(reinterpret_cast<char*>(i_buf.data()), nnz_column * sizeof(uint32_t));
        if (!i_stream_) throw std::runtime_error("Error reading i file: " + iFile());
    } else {
        
        using T = int16_t;
        static constexpr T DELTA_SENTINEL = std::numeric_limits<T>::max();
        // static constexpr T DELTA_SENTINEL = 4;
        std::vector<T> deltas(nnz_column);

        i_stream_.seekg(start * sizeof(T), std::ios::beg);
        i_stream_.read(reinterpret_cast<char*>(deltas.data()), nnz_column * sizeof(T));
        if (!i_stream_) throw std::runtime_error("Error reading delta stream");

        uint64_t delta;
        if (deltas[0] == DELTA_SENTINEL) {
            // overflow → fetch true delta from COO
            delta = I_.get(0, column);
        } else {
            delta = static_cast<uint64_t>(deltas[0]);
        }
        i_buf[0] = static_cast<uint32_t>(column - delta);

        for (size_t idx = 1; idx < nnz_column; ++idx) {
            if (deltas[idx] == DELTA_SENTINEL) {
                delta = I_.get(idx, column);
            } else {
                delta = static_cast<uint64_t>(deltas[idx]);
            }
            i_buf[idx] = i_buf[idx - 1] + static_cast<uint32_t>(delta);
        }
    }

    return i_buf;
}

template <typename XType>
void read_and_scale(std::istream& in, std::vector<float>& vals, size_t nnz, int64_t scale) {
    std::vector<XType> raw(nnz);
    in.read(reinterpret_cast<char*>(raw.data()), nnz * sizeof(XType));
    if (!in) throw std::runtime_error("Error reading input stream");
    constexpr XType NA_SENTINEL = std::numeric_limits<XType>::min();
    for (size_t k = 0; k < nnz; ++k)
        vals[k] = (raw[k] == NA_SENTINEL)
                    ? std::numeric_limits<float>::quiet_NaN()
                    : static_cast<float>(raw[k]) / scale;
}

std::vector<float> LDZipMatrix::get_x(uint32_t column, Stat stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    const auto& p = get_p(); 
    uint64_t nnz = p[column + 1] - p[column];
    std::vector<float> x_buf(nnz);

    std::fstream& x_stream = x_streams_[stat];
    x_stream.clear();    
    if(bits_==Bits::B99)    
    {
        x_stream.seekg(p[column] * sizeof(float), std::ios::beg);
        x_stream.read(reinterpret_cast<char*>(x_buf.data()), nnz * sizeof(float));
    }
    else
    {
        size_t bytes_per_val;
        switch (bits_) {
            case Bits::B8:  bytes_per_val = sizeof(int8_t);  break;
            case Bits::B16: bytes_per_val = sizeof(int16_t); break;
            case Bits::B32: bytes_per_val = sizeof(int32_t); break;
            default: throw std::runtime_error("Unsupported bits value");
        }

        x_stream.seekg(p[column] * bytes_per_val, std::ios::beg);
        int64_t scale = (1LL << (bits_to_int(bits_) - 1)) - 1;
        switch (bits_) {
            case Bits::B8:  read_and_scale<int8_t>(x_stream, x_buf, nnz, scale);  break;
            case Bits::B16: read_and_scale<int16_t>(x_stream, x_buf, nnz, scale); break;
            case Bits::B32: read_and_scale<int32_t>(x_stream, x_buf, nnz, scale); break;
            default: throw std::runtime_error("Unsupported bits value");
        }

    }
    if (!x_stream) throw std::runtime_error("Error reading x file: " + xFile(stat));
    return x_buf;
}


// --- Element Access ---
float LDZipMatrix::getValue(size_t row, size_t col, const Stat& stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    if (row >= nrows_ || col >= ncols_)
        throw std::out_of_range("Index out of bounds: (" + std::to_string(row) + "," + std::to_string(col) + ")");

    auto i_buf = get_i(col);
    auto x_buf = get_x(col, stat);

    for (size_t k = 0; k < i_buf.size(); ++k) {
        if(row == i_buf[k])
            return x_buf[k];
    }

    return 0.0f; 
}


// --- Column/Row Access ---
void LDZipMatrix::getColumn(uint32_t column, std::vector<float>& value, const Stat& stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    value.assign(nrows_, 0.0f);  

    auto i_buf = get_i(column);
    auto x_buf = get_x(column, stat);

    for (size_t k = 0; k < i_buf.size(); ++k) {
        value[i_buf[k]] = x_buf[k];
    }
}

void LDZipMatrix::getColumn(size_t column,
                            size_t row_start,
                            size_t row_end,
                            double* values,
                            const Stat& stat) const { // unsafe, but upto caller to make sure valus is allocated
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    if (row_start > row_end || row_end >= nrows_) {
        throw std::out_of_range("Invalid row range");
    }

    size_t len = row_end - row_start + 1;
    std::fill(values, values + len, 0.0);

    auto i_buf = get_i(column);
    auto x_buf = get_x(column, stat);

    for (size_t k = 0; k < i_buf.size(); ++k) {
        size_t row = i_buf[k];
        if (row >= row_start && row <= row_end) {
            values[row - row_start] = static_cast<double>(x_buf[k]);
        }
    }
}

void LDZipMatrix::getColumn(size_t column,
                            const int* rows, size_t nrows_sel,
                            double* values,
                            const Stat& stat) const { // unsafe, but upto caller to make sure valus is allocated
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    validate_sorted(rows, nrows_sel, "rows");

    std::fill(values, values + nrows_sel, 0.0);

    auto i_buf = get_i(column);
    auto x_buf = get_x(column, stat);

    size_t r_ptr = 0;
    for (size_t k = 0; k < i_buf.size() && r_ptr < nrows_sel; ++k) {
        size_t row = i_buf[k];

        while (r_ptr < nrows_sel && static_cast<size_t>(rows[r_ptr]) < row) {
                ++r_ptr;
            }
        // Handle all duplicate row indices
        while (r_ptr < nrows_sel && static_cast<size_t>(rows[r_ptr]) == row) {
                values[r_ptr] = static_cast<double>(x_buf[k]);
                ++r_ptr;
            }
    }
}

void LDZipMatrix::getColumnRaw(size_t column,
                            std::vector<size_t> rows,
                            std::vector<size_t> &new_i, std::vector<float> &new_x,
                            const Stat& stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    if (column >= ncols_)
        throw std::out_of_range("Invalid column index");

    validate_sorted(rows, "rows");

    auto i_buf = get_i(column);
    auto x_buf = get_x(column, stat);

    size_t r_ptr = 0;
    size_t nr_ptr = 0;
    new_i.resize(rows.size());
    new_x.resize(rows.size());
    for (size_t k = 0; k < i_buf.size() && r_ptr < rows.size(); ++k) {
        size_t row = static_cast<size_t>(i_buf[k]);

        while (r_ptr < rows.size() && static_cast<size_t>(rows[r_ptr]) < row) {
            ++r_ptr;
        }

        // Handle all duplicate row indices
        while (r_ptr < rows.size() && static_cast<size_t>(rows[r_ptr]) == row) {
            new_i[nr_ptr] = r_ptr;
            new_x[nr_ptr++] = x_buf[k];
            ++r_ptr;
        }
    }
    new_i.resize(nr_ptr);
    new_x.resize(nr_ptr);

}

// --- Sub Matrix Access
void LDZipMatrix::getSubMatrix(size_t col_start,
                               size_t col_end,
                               size_t row_start,
                               size_t row_end,
                               double* values,
                               const Stat& stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");
    if (col_start > col_end || col_end >= ncols_)
        throw std::out_of_range("Invalid column range");
    if (row_start > row_end || row_end >= nrows_)
        throw std::out_of_range("Invalid row range");

    size_t nrows_sel = row_end - row_start + 1;
    size_t ncols_sel = col_end - col_start + 1;

    // Fill with zeros first
    std::fill(values, values + (nrows_sel * ncols_sel), 0.0);

    for (uint32_t col = col_start; col <= col_end; ++col) {
        auto i_buf = get_i(col);
        auto x_buf = get_x(col, stat);

        size_t col_idx = col - col_start;
        for (size_t k = 0; k < i_buf.size(); ++k) {
            uint32_t row = i_buf[k];
            if (row >= row_start && row <= row_end) {
                values[(row - row_start) + col_idx * nrows_sel] =
                    static_cast<double>(x_buf[k]);
            }
        }
    }
}

void LDZipMatrix::getSubMatrix(const int* cols, size_t ncols_sel,
                               const int* rows, size_t nrows_sel,
                               double* values, const Stat& stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    validate_sorted(rows, nrows_sel, "rows");

    std::fill(values, values + (nrows_sel * ncols_sel), 0.0);

    for (size_t j = 0; j < ncols_sel; ++j) {
        size_t col = static_cast<size_t>(cols[j]); // adjust if R is 1-based
        if (col >= ncols_)
            throw std::out_of_range("Invalid column index");

        auto i_buf = get_i(col);
        auto x_buf = get_x(col, stat);

        size_t r_ptr = 0;
        for (size_t k = 0; k < i_buf.size() && r_ptr < nrows_sel; ++k) {
            size_t row = static_cast<size_t>(i_buf[k]);

            while (r_ptr < nrows_sel && static_cast<size_t>(rows[r_ptr]) < row) {
                ++r_ptr;
            }

            // Handle all duplicate row indices
            while (r_ptr < nrows_sel && static_cast<size_t>(rows[r_ptr]) == row) {
                values[r_ptr + j * nrows_sel] = static_cast<double>(x_buf[k]);
                ++r_ptr;
            }
        }
    }
}

// --- Pairwise Access

void LDZipMatrix::getPairwise(const int* cols,
                            const int* rows,
                            size_t npairs,
                            double* values,
                            const Stat& stat) const {

    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");

    // If sorted by column, this will optimize by grouping same column queries
    // If not sorted, it still works correctly (just less efficiently)
    size_t start = 0;
    while (start < npairs) {

        size_t col = static_cast<size_t>(cols[start]);
        size_t end = start + 1;

        // group same column queries
        while (end < npairs && cols[end] == cols[start])
            ++end;

        size_t block_size = end - start;

        getColumn(col, rows + start, block_size, values + start, stat );

        start = end;
    }
}

// --- Get Neighbors ---
std::vector<uint32_t> LDZipMatrix::get_neighbors(uint32_t column, double abs_threshold, const Stat& stat) const {
    if (!has_stat_[stat])
        throw std::invalid_argument(" Stats Type [" + stat_to_string(stat) + "] not available");
    const auto& i_buf = get_i(column);
    const auto& x_buf = get_x(column, stat);
    std::vector<uint32_t> neighbors;
    neighbors.reserve(i_buf.size()); 
    for (size_t k = 0; k < i_buf.size(); ++k) {
        size_t row = i_buf[k];
        double val = static_cast<double>(x_buf[k]);
        if (std::abs(val) >= abs_threshold) {
            neighbors.push_back(row);
        }
    }
    return neighbors;
} 


// --- Variant Metadata
void LDZipMatrix::readVariants(const std::string& snp_file) {
    variants_.clear();

    std::ifstream vars_in(snp_file);
    if (!vars_in) {
        throw std::runtime_error("Cannot open SNP vars file : " + snp_file);
    }

    std::string line;
    while (std::getline(vars_in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }
        if (fields.size() < 5)
            throw std::runtime_error("SNP file has fewer than 5 columns: " + line);

        variants_.push_back(Variant(fields[2], fields[3], fields[4]));
    }

    if (variants_.size() != nrows_) {
        throw std::runtime_error("Variant count does not match matrix row count");
    }
}

MetaInfo LDZipMatrix::metaInfo() const{

    MetaInfo meta(nrows_, ncols_, nnz_, bits(), format_, version_);
    for (Stat s : stats_available_) meta.has_stat[s] = true;
    return meta;
}

// --- Decompress ---
void LDZipMatrix::writeColumnBinary(     std::ofstream& out, 
                                                const std::vector<float>& value) {
    out.write(reinterpret_cast<const char*>(value.data()), nrows_ * sizeof(float));
    if (!out) {
        throw std::runtime_error("Error writing column to binary file");
    }
}

void LDZipMatrix::writeColumnTabular(
        std::ofstream& out, 
        size_t column){
    if (variants_.empty()) {
        throw std::runtime_error(" variants have not been initialized. contact author with command line ");
    }

    auto i_buf = get_i(column);
    auto x_bufs = EnumArray<std::vector<float>, Stat>();

    for (Stat s : stats_available_)
        x_bufs[s] = get_x(column, s);

    for (size_t idx = 0; idx < i_buf.size(); ++idx) {
        auto row = i_buf[idx];
        if (row <= column) continue;

        const Variant& va = variants_[column];
        const Variant& vb = variants_[row];

        out << va.id()  << '\t' << va.ref()  << '\t' << va.alt()  << '\t'
            << vb.id()  << '\t' << vb.ref()  << '\t' << vb.alt();

        for (Stat s : stats_available_) 
            out << '\t' << x_bufs[s][idx];
        out << '\n';
    }

    if (!out) {
        throw std::runtime_error("Error writing column to tabular file");
    }
}

void LDZipMatrix::decompress(const std::string& out_prefix,
                                    const std::string& type) {
    std::vector<float> x_buf(nrows_);

    if (type == "binary") {
        for (Stat s : stats_available_) {
                std::ofstream out(out_prefix + "." + stat_to_string(s) + ".bin", std::ios::binary);
                if (!out) {
                    throw std::runtime_error("Error opening binary output file: " + out_prefix);
                }
                for (size_t col = 0; col < ncols_; ++col) {
                    getColumn(col, x_buf, s);
                    writeColumnBinary(out, x_buf);
                }
            }    
    }
    else if (type == "tabular") {
        std::ofstream out(out_prefix + ".vcor");
        if (!out) {
            throw std::runtime_error("Error opening tabular output file: " + out_prefix + ".vcor");
        }
        out << "#ID_A\tREF_A\tALT_A\tID_B\tREF_B\tALT_B";
        for (Stat s : stats_available_) 
            out << "\t" << stat_to_string(s);
        out << "\n";
        for (size_t col = 0; col < ncols_; ++col) {
            writeColumnTabular(out, col);
        }
    }
    else {
        throw std::runtime_error("Unsupported decompress type: " + type);
    }
}


}
