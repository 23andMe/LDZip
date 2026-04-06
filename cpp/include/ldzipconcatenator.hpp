#pragma once
#include "ldzipcompressor.hpp"
#include "ldzipmatrix.hpp"
#include <string>
#include <vector>

namespace ldzip {

// Per-chunk boundary information for overlapping concat.
// Local indices are 0-based positions within the chunk's own variant list.
struct ChunkBoundary {
    size_t n_variants;            // total variants in this chunk
    size_t first_overlap_end;     // local excl-end of overlap with chunk[i-1]  (0 if first chunk)
    size_t second_overlap_start;  // local start of overlap with chunk[i+1]     (n_variants if last chunk)
    size_t global_start;          // index in merged array where this chunk's window starts
    size_t global_end;            // index in merged array where this chunk's window ends (exclusive)
};

struct OverlapVariantInfo {
    size_t total_variants;
    std::vector<ChunkBoundary> chunks;
};

// Builds variant boundary information for concatenation
// If overlapping=true: Reads variant files and validates overlaps match exactly
// If overlapping=false: Constructs non-overlapping boundaries (treats chunks as independent)
OverlapVariantInfo read_overlapping_variant_order(const std::vector<std::string>& var_files, bool overlapping);

class LDZipConcatenator {
public:
    LDZipConcatenator(size_t nrows, size_t ncols, MatrixFormat format, const std::vector<Stat>& stats, Bits bits, const std::string& prefix);
    void process_overlap_columns(const LDZipMatrix& current_mat, const LDZipMatrix& next_mat, const ChunkBoundary& current_chunk, const ChunkBoundary& next_chunk);
    void process_exclusive_columns(const LDZipMatrix& current_mat, const ChunkBoundary& current_chunk);
    void close();

private:
    // Compressor used for merging overlapping regions
    LDZipCompressor overlap_merger_;

    // Helper methods
    void write_initial_p_zero();  // Write initial 0 to p_stream_
    void init_p_for_column(uint32_t col, uint64_t value);
    uint64_t get_p_value(uint32_t col) const;
    void push_decoded_column(uint32_t cidx);
    void writeMergedColumn();

    // Helper functions for bulk file copying
    template <typename T>
    void copy_binary_file_with_offset(const std::string& in_file, std::fstream& out_stream, std::streamoff read_pos, T add_value, size_t byte_size, size_t n, size_t chunk_size = 1000000);

    template <typename T>
    void copy_binary_file_with_offset_file(const std::string& in_file, const std::string& out_file, std::streamoff read_pos, T add_value, size_t byte_size, size_t n, size_t chunk_size = 1000000);

    std::pair<uint64_t, uint64_t> get_I_range(const std::string& f, size_t a, size_t b);

    // Helper to clear merge buffers
    void clear_merge_buffers(const LDZipMatrix& mat);

    // Helper to extract filtered column data into above/below buffers
    void extract_above_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t threshold, uint32_t global_offset);
    void extract_below_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t threshold, uint32_t global_offset);

    // Helper to merge above and below buffers into merged buffers
    void merge_above_below_buffers(const LDZipMatrix& mat);

    // Overlap merge buffers - reused across columns to avoid repeated allocation
    std::vector<uint32_t> above_rows_, above_idx_, below_rows_, below_idx_, merged_rows_;
    EnumArray<std::vector<float>, Stat> above_x_, below_x_, merged_x_;

    // Track total nnz internally as we process chunks
    uint64_t current_nnz_ = 0;
};

} // namespace ldzip
