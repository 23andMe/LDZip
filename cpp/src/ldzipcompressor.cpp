#include "ldzipcompressor.hpp"

namespace ldzip
{

LDZipCompressor::LDZipCompressor(size_t nrows,
                         size_t ncols,
                         MatrixFormat format,
                         const std::vector<Stat>& stats,
                         Bits bits,
                         const std::string& prefix,
                         Mode mode)
    :   m_(nrows, ncols, format, stats, bits, prefix),
        mode_(mode),
        diag_vals{1.0f} {
    
    p_stream_.open(m_.pFile(), std::ios::out | std::ios::binary | std::ios::trunc);
    i_stream_.open(m_.iFile(), std::ios::out | std::ios::binary | std::ios::trunc);
    for (Stat s : m_.stats_available_) 
        x_streams_[s].open(m_.xFile(s), std::ios::out | std::ios::binary | std::ios::trunc);

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
                         Mode mode)
    : LDZipCompressor(nrows, ncols, format, std::vector<Stat>{stat}, bits, prefix, mode) {    
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
        
        if(m_.version_ == "1.1")
        {
            i_stream_.write(
                reinterpret_cast<const char*>(i_col.data()), static_cast<std::streamsize>(i_col.size() * sizeof(uint32_t))
            );
        } else
        {
            
            using T = int16_t;
            static constexpr T DELTA_SENTINEL = std::numeric_limits<T>::max();
            // static constexpr T DELTA_SENTINEL = 4;

            std::vector<T> deltas(i_col.size());
            uint64_t delta = static_cast<int64_t>(active_column_) - static_cast<int64_t>(i_col[0]);
            if (delta >= DELTA_SENTINEL ) {
                m_.I_.push(static_cast<uint32_t>(delta));
                deltas[0] = DELTA_SENTINEL;
            }else {
                deltas[0] = static_cast<T>(delta);
            }

            for (size_t idx = 1; idx < i_col.size(); ++idx) {
                delta = static_cast<int64_t>(i_col[idx]) - static_cast<int64_t>(i_col[idx - 1]);
                if (delta >= DELTA_SENTINEL) {
                    m_.I_.push(static_cast<uint32_t>(delta));
                    deltas[idx] = DELTA_SENTINEL;
                }else {
                    deltas[idx] = static_cast<T>(delta);
                }
            }

            i_stream_.write(reinterpret_cast<const char*>(deltas.data()),
                            static_cast<std::streamsize>(deltas.size() * sizeof(T)));

        }
    }
    m_.I_.flush_column();
    i_col.clear();
    i_col.shrink_to_fit();
}

void LDZipCompressor::write_x(Stat& s){

    auto& x_col = m_.xs_[s][active_column_];
    std::fstream& x_out = x_streams_[s];
    int64_t scale = (1LL << (m_.bits() - 1)) - 1;

    switch (m_.bits_) {
        case Bits::B8:
            write_scaled_buffer<int8_t>(x_out, x_col, scale);
            break;
        case Bits::B16:
            write_scaled_buffer<int16_t>(x_out, x_col, scale);
            break;
        case Bits::B32:
            write_scaled_buffer<int32_t>(x_out, x_col, scale);
            break;
        case Bits::B99:
            write_scaled_buffer<float>(x_out, x_col, scale);
            break;
        default:
            throw std::runtime_error("Unsupported bits value");
    }

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
    for (Stat s : m_.stats_available_) 
        x_streams_[s].close();
    i_stream_.close();
    p_stream_.close();

    write_metadata_json(m_.metaFile(), m_.metaInfo());

}


} // namespace ldzip


  


  