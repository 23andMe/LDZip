#include "ldzipcompressor.hpp"

namespace ldzip
{

LDZipCompressor::LDZipCompressor(size_t nrows,
                         size_t ncols,
                         MatrixFormat format,
                         const std::vector<Stat>& stats,
                         Bits bits,
                         const std::string& prefix,
                         Mode mode,
                         size_t chunk_size)
    :   m_(nrows, ncols, format, stats, bits, prefix, chunk_size),
        mode_(mode),
        diag_vals{1.0f} {

    if (chunk_size == 0)
        throw std::runtime_error("chunk_size must be > 0 for v3.0");

    p_stream_.open(m_.pFile(), std::ios::out | std::ios::binary | std::ios::trunc);

    i_chunked_writer_ = std::make_unique<ChunkedWriter>(m_.iFile(), m_.iIndexFile(), chunk_size);
    for (Stat s : m_.stats_available_)
        x_chunked_writers_[s] = std::make_unique<ChunkedWriter>(m_.xFile(s), m_.xIndexFile(s), chunk_size);

    diag_vals[Stat::D] = std::numeric_limits<float>::quiet_NaN();
    active_column_ = -1;
    if(mode == Mode::ColumnStream)
    {
        x_buffer.resize(nrows);
        active_column_ = 0;
    }
}

LDZipCompressor::LDZipCompressor(size_t nrows,
                         size_t ncols,
                         MatrixFormat format,
                         Stat stat,
                         Bits bits,
                         const std::string& prefix,
                         Mode mode,
                         size_t chunk_size)
    : LDZipCompressor(nrows, ncols, format, std::vector<Stat>{stat}, bits, prefix, mode, chunk_size) {
}

void LDZipCompressor::push_value(
                uint32_t ridx,
                uint32_t cidx,
                const EnumArray<float, Stat>& values){
    
    if (static_cast<int>(cidx) < active_column_)
        throw std::runtime_error("vcor table NOT column by column, each column cut off at the diagonal.");
    if (cidx >= ridx)
        throw std::runtime_error("vcor table NOT column by column, each column cut off at the diagonal.");
    
    if (static_cast<int>(cidx) != active_column_){

        if( active_column_ >= 0 ) 
            writeActiveColumn();
        active_column_++;

        while( active_column_ < static_cast<int>(cidx) ){
            addTrivialColumn();
            writeActiveColumn();
            active_column_++;
        }
        push_value_(cidx, cidx, diag_vals);
    
    }
    push_value_(ridx, cidx, values);
    push_value_(cidx, ridx, values);
}

void LDZipCompressor::push_column_raw(
                        uint32_t cidx,
                        const std::vector<float>& values,
                        const std::vector<size_t>& rindices,
                        Stat stat) {

    EnumArray<float, Stat> Statvalues(-999.0f);

    size_t start_r = 0;
    size_t end_r = rindices.size() - 1;
    if (m_.format_ == MatrixFormat::UPPER) { start_r = cidx; }

    for (size_t r = start_r; r <= end_r; ++r){
        size_t row = rindices[r];
        float val = values[r];
        if (std::isnan(val)) {
            throw std::runtime_error("NaN encountered in values despite prior filtering. Genotype data must have missing values");
        }
        Statvalues[stat] = val;
        push_value_(row, cidx, Statvalues);
    }

    writeActiveColumn();
    active_column_++;

}

void LDZipCompressor::push_column(
                        uint32_t cidx,
                        const std::vector<float>& values,
                        const std::vector<size_t>& keep,
                        float min,
                        Stat stat) {

    EnumArray<float, Stat> Statvalues(-999.0f);

    size_t start_r = 0;
    size_t end_r = m_.nrows_ - 1;
    if (m_.format_ == MatrixFormat::UPPER) { start_r = cidx; }

    for (size_t r = start_r; r <= end_r; ++r){
        size_t row = keep[r];
        float val = values[row];
        if (!std::isnan(val) && std::abs(val) < min) continue;

        Statvalues[stat] = val;
        push_value_(r, cidx, Statvalues);
    }

    writeActiveColumn();
    active_column_++;

}

// Encode x values to byte buffer for chunked writing
template <typename T>
std::vector<uint8_t> encode_scaled_buffer(const std::vector<float>& x_col, int64_t scale) {
    std::vector<T> buf(x_col.size());
    for (size_t i = 0; i < x_col.size(); ++i) {
        float v = x_col[i];
        if (std::is_same<T, float>::value) {
            buf[i] = v;
        } else {
            if (std::isnan(v))
                buf[i] = std::numeric_limits<T>::min();
            else{
                double prod = static_cast<double>(v) * static_cast<double>(scale);
                buf[i] = static_cast<T>(std::llround(prod));
            }
        }
    }
    std::vector<uint8_t> result(buf.size() * sizeof(T));
    std::memcpy(result.data(), buf.data(), result.size());
    return result;
}

template <typename T, typename Stream>
void write_scaled_buffer(Stream& x_out, const std::vector<float>& x_col, int64_t scale) {
    std::vector<T> buf(x_col.size());
    for (size_t i = 0; i < x_col.size(); ++i) {
        float v = std::clamp(x_col[i], -1.0f, 1.0f);
        if constexpr (std::is_same_v<T, float>) {
            // B99: raw float (fine to add NaN values here)
            buf[i] = v;
        } else {
            if (std::isnan(v))
                buf[i] = std::numeric_limits<T>::min();
            else{
                double prod = static_cast<double>(v) * static_cast<double>(scale);
                buf[i] = static_cast<T>(std::llround(prod));
            }
        }
    }
    x_out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size() * sizeof(T)));
}

void LDZipCompressor::writeActiveColumn(){

    // Update p_vector before writing out column
    m_.p_[active_column_ + 1] = m_.p_[active_column_] + m_.i_[active_column_].size();

    for (Stat s : All_Stats()) if (m_.has_stat_[s]) write_x(s);
    write_i();
}

void LDZipCompressor::write_i() {
    std::vector<uint32_t>& i_col = m_.i_[active_column_];

    if (!i_col.empty()) {
        // v3.0: Encode as int32_t deltas (first element from diagonal, rest as successive differences)
        std::vector<int32_t> deltas(i_col.size());
        deltas[0] = static_cast<int32_t>(i_col[0]) - static_cast<int32_t>(active_column_);
        for (size_t idx = 1; idx < i_col.size(); ++idx) {
            deltas[idx] = static_cast<int32_t>(i_col[idx]) - static_cast<int32_t>(i_col[idx - 1]);
        }
        i_chunked_writer_->appendColumn(deltas.data(), deltas.size() * sizeof(int32_t));
    }

    i_col.clear();
    i_col.shrink_to_fit();
}

void LDZipCompressor::write_x(Stat& s){

    auto& x_col = m_.xs_[s][active_column_];
    int64_t scale = (1LL << (m_.bits() - 1)) - 1;

    std::vector<uint8_t> encoded_data;
    switch (m_.bits_) {
        case Bits::B8:
            encoded_data = encode_scaled_buffer<int8_t>(x_col, scale);
            break;
        case Bits::B16:
            encoded_data = encode_scaled_buffer<int16_t>(x_col, scale);
            break;
        case Bits::B32:
            encoded_data = encode_scaled_buffer<int32_t>(x_col, scale);
            break;
        case Bits::B99:
            encoded_data = encode_scaled_buffer<float>(x_col, scale);
            break;
        default:
            throw std::runtime_error("Unsupported bits value");
    }

    if (!encoded_data.empty())
        x_chunked_writers_[s]->appendColumn(encoded_data.data(), encoded_data.size());

    x_col.clear();
    x_col.shrink_to_fit();
}

void LDZipCompressor::write_p(){
    if (!m_.p_.empty()) {
        p_stream_.write(reinterpret_cast<const char*>(m_.p_.data()),
                        static_cast<std::streamsize>(m_.p_.size() * sizeof(uint64_t)));
    }
}

void LDZipCompressor::addTrivialColumn(){

    for (Stat s : m_.stats_available_) {
        auto& x = m_.xs_[s][active_column_];
        if (s == Stat::D)
            x.push_back(std::numeric_limits<float>::quiet_NaN());
        else
            x.push_back(1.0f);
    }

    m_.i_[active_column_].push_back(active_column_);
    m_.nnz_++;
}

void LDZipCompressor::push_value_(
                uint32_t ridx,
                uint32_t cidx,
                const EnumArray<float, Stat>& values){
                    
    for (Stat s : m_.stats_available_) {
        float value = values[s];
        auto& x = m_.xs_[s][cidx];
        x.push_back(value);
    }

    m_.i_[cidx].push_back(ridx);
    m_.nnz_++;
}

void LDZipCompressor::stream_close()
{
    if(mode_ == Mode::ValueStream)
    {
        if( active_column_ >= 0 )
            writeActiveColumn();
        active_column_++;

        while(active_column_ < static_cast<int>(m_.ncols_))
        {
            addTrivialColumn();
            writeActiveColumn();
            active_column_++;
        }
    }

    write_p();
    p_stream_.close();

    // Close all chunked writers (v3.0)
    i_chunked_writer_->close();
    for (Stat s : m_.stats_available_) {
        x_chunked_writers_[s]->close();
    }

    write_metadata_json(m_.metaFile(), m_.metaInfo());

}


} // namespace ldzip


  


  