#pragma once
#include <zstd.h>
#include <vector>
#include <fstream>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <memory>
#include <unordered_map>
#include <algorithm>

namespace ldzip {

    // Chunked writer for zstd-compressed column data
    class ChunkedWriter {
        private:
            std::ofstream chunks_file_;
            std::ofstream index_file_;
            std::vector<uint8_t> buffer_;
            std::vector<uint8_t> compressed_buffer_;  // Reusable compression buffer
            size_t chunk_size_;          // Number of columns per chunk
            size_t current_chunk_cols_;  // Columns buffered in current chunk
            uint64_t current_offset_;    // Current byte offset in chunks file
            size_t total_columns_;       // Total columns written
            std::vector<uint64_t> byte_offsets_;   // Byte offsets for each chunk boundary
            std::vector<uint64_t> column_counts_;  // Cumulative column counts at each boundary
            int compression_level_;
            bool closed_ = false;

        public:
            ChunkedWriter(const std::string& chunks_path,
                         const std::string& index_path,
                         size_t chunk_size,
                         int compression_level = 6)
                : chunk_size_(chunk_size),
                  current_chunk_cols_(0),
                  current_offset_(0),
                  total_columns_(0),
                  compression_level_(compression_level) {

                chunks_file_.open(chunks_path, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!chunks_file_) {
                    throw std::runtime_error("Cannot open chunks file for write: " + chunks_path);
                }

                index_file_.open(index_path, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!index_file_) {
                    throw std::runtime_error("Cannot open index file for write: " + index_path);
                }

                // Pre-allocate arrays (estimate ~200 chunks for typical use)
                byte_offsets_.reserve(200);
                column_counts_.reserve(200);
                byte_offsets_.push_back(0);
                column_counts_.push_back(0);
            }

            // Append raw bytes for one column to the buffer
            void appendColumn(const void* data, size_t size) {
                const uint8_t* bytes = static_cast<const uint8_t*>(data);
                buffer_.insert(buffer_.end(), bytes, bytes + size);
                current_chunk_cols_++;
                total_columns_++;

                // Flush when chunk is full
                if (current_chunk_cols_ >= chunk_size_) {
                    flush();
                }
            }

            // Compress and write current buffer
            void flush() {
                if (buffer_.empty()) return;

                // Compress buffer with zstd (reuse compression buffer)
                size_t compressed_bound = ZSTD_compressBound(buffer_.size());
                if (compressed_buffer_.size() < compressed_bound) {
                    compressed_buffer_.resize(compressed_bound);
                }

                size_t compressed_size = ZSTD_compress(
                    compressed_buffer_.data(), compressed_buffer_.size(),
                    buffer_.data(), buffer_.size(),
                    compression_level_
                );

                if (ZSTD_isError(compressed_size)) {
                    throw std::runtime_error("ZSTD compression failed: " +
                                           std::string(ZSTD_getErrorName(compressed_size)));
                }

                // Write compressed data
                chunks_file_.write(reinterpret_cast<const char*>(compressed_buffer_.data()), compressed_size);
                current_offset_ += compressed_size;

                // Track boundary in arrays (will write to index on close)
                byte_offsets_.push_back(current_offset_);
                column_counts_.push_back(total_columns_);

                // Clear buffer
                buffer_.clear();
                current_chunk_cols_ = 0;
            }

            void close() {
                if (closed_) return;
                flush(); // Write any remaining data

                // Compress index: [byte_offsets][column_counts] with zstd
                size_t n_boundaries = byte_offsets_.size();
                size_t index_raw_size = 2 * n_boundaries * sizeof(uint64_t);

                // Reuse compressed_buffer_ for index compression
                size_t compressed_bound = ZSTD_compressBound(index_raw_size);
                if (compressed_buffer_.size() < compressed_bound) {
                    compressed_buffer_.resize(compressed_bound);
                }

                // Write directly into buffer without intermediate vector
                uint64_t* index_ptr = reinterpret_cast<uint64_t*>(buffer_.data());
                if (buffer_.size() < index_raw_size) {
                    buffer_.resize(index_raw_size);
                    index_ptr = reinterpret_cast<uint64_t*>(buffer_.data());
                }

                std::memcpy(index_ptr, byte_offsets_.data(), n_boundaries * sizeof(uint64_t));
                std::memcpy(index_ptr + n_boundaries, column_counts_.data(), n_boundaries * sizeof(uint64_t));

                size_t compressed_size = ZSTD_compress(
                    compressed_buffer_.data(), compressed_bound,
                    buffer_.data(), index_raw_size,
                    compression_level_
                );

                if (ZSTD_isError(compressed_size)) {
                    throw std::runtime_error("Index compression failed: " +
                                           std::string(ZSTD_getErrorName(compressed_size)));
                }

                // Write compressed index
                index_file_.write(reinterpret_cast<const char*>(compressed_buffer_.data()), compressed_size);

                if (chunks_file_.is_open()) chunks_file_.close();
                if (index_file_.is_open()) index_file_.close();
                closed_ = true;
            }

            ~ChunkedWriter() {
                close();
            }
    };

    // Chunked reader with LRU cache for decompressed chunks
    class ChunkedReader {
        private:
            std::string chunks_path_;
            std::string index_path_;
            mutable std::ifstream chunks_file_;
            std::vector<uint64_t> byte_offsets_;   // Loaded from index
            std::vector<uint64_t> column_counts_;  // Loaded from index
            size_t num_chunks_;

            // Reusable buffers to avoid reallocation on every read
            mutable std::vector<uint8_t> compressed_buffer_;
            mutable std::vector<uint8_t> decompressed_buffer_;

            // Cache: chunk_id -> decompressed data
            struct CacheEntry {
                std::vector<uint8_t> data;
                size_t access_count;
            };
            mutable std::unordered_map<size_t, CacheEntry> cache_;
            mutable size_t cache_access_counter_;
            static constexpr size_t MAX_CACHE_SIZE = 10; // Cache up to 10 chunks

        public:
            ChunkedReader(const std::string& chunks_path,
                         const std::string& index_path,
                         size_t /* chunk_size unused */)
                : chunks_path_(chunks_path),
                  index_path_(index_path),
                  cache_access_counter_(0) {

                // Read and decompress index file
                std::ifstream index_file(index_path_, std::ios::binary | std::ios::ate);
                if (!index_file) {
                    throw std::runtime_error("Cannot open index file: " + index_path_);
                }

                size_t compressed_size = index_file.tellg();
                index_file.seekg(0);

                // Use compressed_buffer_ for reading (reuse member buffer)
                if (compressed_buffer_.size() < compressed_size) {
                    compressed_buffer_.resize(compressed_size);
                }
                index_file.read(reinterpret_cast<char*>(compressed_buffer_.data()), compressed_size);
                index_file.close();

                // Decompress index
                unsigned long long decompressed_size = ZSTD_getFrameContentSize(
                    compressed_buffer_.data(), compressed_size
                );

                if (decompressed_size == ZSTD_CONTENTSIZE_ERROR ||
                    decompressed_size == ZSTD_CONTENTSIZE_UNKNOWN) {
                    throw std::runtime_error("Cannot determine index decompressed size");
                }

                // Use decompressed_buffer_ for decompression (reuse member buffer)
                if (decompressed_buffer_.size() < decompressed_size) {
                    decompressed_buffer_.resize(decompressed_size);
                }

                size_t result = ZSTD_decompress(
                    decompressed_buffer_.data(), decompressed_size,
                    compressed_buffer_.data(), compressed_size
                );

                if (ZSTD_isError(result)) {
                    throw std::runtime_error("Index decompression failed: " +
                                           std::string(ZSTD_getErrorName(result)));
                }

                // Split into byte_offsets and column_counts
                size_t n_entries = decompressed_size / sizeof(uint64_t);
                size_t n_boundaries = n_entries / 2;

                uint64_t* data_ptr = reinterpret_cast<uint64_t*>(decompressed_buffer_.data());
                byte_offsets_.assign(data_ptr, data_ptr + n_boundaries);
                column_counts_.assign(data_ptr + n_boundaries, data_ptr + n_entries);

                num_chunks_ = n_boundaries - 1;

                chunks_file_.open(chunks_path_, std::ios::in | std::ios::binary);
                if (!chunks_file_) {
                    throw std::runtime_error("Cannot open chunks file: " + chunks_path_);
                }
            }

            // Read and decompress a specific chunk
            const std::vector<uint8_t>& readChunk(size_t chunk_id) const {
                // Check cache first
                auto it = cache_.find(chunk_id);
                if (it != cache_.end()) {
                    it->second.access_count = ++cache_access_counter_;
                    return it->second.data;
                }

                // Get chunk byte range from loaded index
                uint64_t start_offset = byte_offsets_[chunk_id];
                uint64_t end_offset = byte_offsets_[chunk_id + 1];
                size_t compressed_size = end_offset - start_offset;

                // Read compressed data (reuse buffer to avoid reallocation)
                if (compressed_buffer_.size() < compressed_size) {
                    compressed_buffer_.resize(compressed_size);
                }
                chunks_file_.clear();
                chunks_file_.seekg(start_offset);
                chunks_file_.read(reinterpret_cast<char*>(compressed_buffer_.data()), compressed_size);

                // Decompress
                unsigned long long decompressed_size = ZSTD_getFrameContentSize(
                    compressed_buffer_.data(), compressed_size
                );

                if (decompressed_size == ZSTD_CONTENTSIZE_ERROR ||
                    decompressed_size == ZSTD_CONTENTSIZE_UNKNOWN) {
                    throw std::runtime_error("Cannot determine decompressed size");
                }

                // Reuse decompression buffer to avoid reallocation
                if (decompressed_buffer_.size() < decompressed_size) {
                    decompressed_buffer_.resize(decompressed_size);
                }
                size_t result = ZSTD_decompress(
                    decompressed_buffer_.data(), decompressed_size,
                    compressed_buffer_.data(), compressed_size
                );

                if (ZSTD_isError(result)) {
                    throw std::runtime_error("ZSTD decompression failed: " +
                                           std::string(ZSTD_getErrorName(result)));
                }

                // Add to cache (evict oldest if full)
                if (cache_.size() >= MAX_CACHE_SIZE) {
                    // Evict LRU entry
                    auto oldest = cache_.begin();
                    for (auto it = cache_.begin(); it != cache_.end(); ++it) {
                        if (it->second.access_count < oldest->second.access_count) {
                            oldest = it;
                        }
                    }
                    cache_.erase(oldest);
                }

                // Copy only the used portion of the buffer into cache
                std::vector<uint8_t> cached_data(decompressed_buffer_.begin(),
                                                  decompressed_buffer_.begin() + decompressed_size);
                cache_[chunk_id] = {std::move(cached_data), ++cache_access_counter_};
                return cache_[chunk_id].data;
            }

            // Find which chunk contains the given column
            size_t getChunkForColumn(size_t col) const {
                // Binary search in column_counts_
                auto it = std::upper_bound(column_counts_.begin(), column_counts_.end(), col);
                if (it == column_counts_.begin()) {
                    throw std::runtime_error("Column " + std::to_string(col) + " out of range");
                }
                return std::distance(column_counts_.begin(), it) - 1;
            }

            // Get the starting column index of a chunk
            size_t getChunkStartColumn(size_t chunk_id) const {
                return column_counts_[chunk_id];
            }

            size_t getNumChunks() const { return num_chunks_; }
            size_t getTotalColumns() const { return column_counts_.back(); }
            uint64_t getByteOffset(size_t chunk_id) const { return byte_offsets_[chunk_id]; }
            uint64_t getChunkSize(size_t chunk_id) const { return byte_offsets_[chunk_id + 1] - byte_offsets_[chunk_id]; }

            void preloadColumns(uint32_t start_col, uint32_t end_col) const {
                size_t start_chunk = getChunkForColumn(start_col);
                size_t end_chunk = getChunkForColumn(end_col);
                uint64_t start_byte = byte_offsets_[start_chunk];
                uint64_t end_byte = byte_offsets_[end_chunk + 1];

                std::vector<char> buffer(64 * 1024 * 1024);
                std::ifstream file(chunks_path_, std::ios::binary);
                if (!file) throw std::runtime_error("Cannot open: " + chunks_path_);
                file.seekg(start_byte);
                uint64_t remaining = end_byte - start_byte;
                while (remaining > 0 && file) {
                    size_t to_read = std::min(remaining, (uint64_t)buffer.size());
                    file.read(buffer.data(), to_read);
                    remaining -= file.gcount();
                }
            }

            ~ChunkedReader() {
                if (chunks_file_.is_open()) chunks_file_.close();
            }
    };

} // namespace ldzip
