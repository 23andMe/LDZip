#pragma once
#include "ldzipmatrix.hpp"
#include <cstdint>
#include <vector>
#include <string>

namespace ldzip {

class LDZipCompressor {
    friend class LDZipConcatenator;
public:
    enum class Mode : uint8_t {
        ValueStream,
        ColumnStream,
        RowStream
    };

    explicit LDZipCompressor    (size_t nrows, size_t ncols, MatrixFormat format, const std::vector<Stat>& stats, Bits bits, const std::string& prefix, Mode mode, size_t chunk_size = 0);
    explicit LDZipCompressor    (size_t nrows, size_t ncols, MatrixFormat format, Stat stats, Bits bits, const std::string& prefix, Mode mode, size_t chunk_size = 0);

    // explicit LDZipCompressor    (LDMatrix& m, Mode mode);

    Mode mode() const;


    // Public writers for compressed file
    void stream_close();
    void push_value(uint32_t ridx, uint32_t cidx, const EnumArray<float, Stat>& values);
    void push_column(uint32_t cidx, const std::vector<float>& values, const std::vector<size_t>& keep, float min, Stat stat);
    void push_column_raw(uint32_t cidx, const std::vector<float>& values, const std::vector<size_t>& indices, Stat stat);


private:

    // --- Private matrix
    LDZipMatrix m_;

    //  --- Private file handlers
    Mode mode_;
    mutable std::fstream p_stream_;

    // v3.0 chunked compression (always used for writing)
    std::unique_ptr<ChunkedWriter> i_chunked_writer_;
    EnumArray<std::unique_ptr<ChunkedWriter>, Stat> x_chunked_writers_;

    // Legacy streams (for concatenator friend access only - not used in v3.0)
    mutable std::fstream i_stream_;
    EnumArray<std::fstream, Stat> x_streams_;

    // Private compressor helpers
    void push_value_(uint32_t ridx, uint32_t cidx, const EnumArray<float, Stat>& values);
    void writeActiveColumn();
    void addTrivialColumn();
    void write_i();
    void write_x(Stat &s);
    void write_p();
    
    EnumArray<float, Stat> diag_vals;  
    int active_column_;
    std::vector<float> x_buffer;

public: 
    const LDZipMatrix& m() const {return m_;}
    LDZipMatrix& m() {return m_;}

};

} // namespace ldzip
