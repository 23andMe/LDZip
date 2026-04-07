#include <fstream>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <string>
#include <algorithm>

namespace ldzip
{
class COO {
    private:
        mutable std::fstream fstream_;
        std::ofstream index_stream_;
        std::string value_file_;   // .io.bin
        std::string index_file_;   // .io.index
        char mode_;

        // Write mode state
        uint64_t index_;

        // Read mode state
        mutable std::vector<uint32_t> column_cache_;  // Current column's values (reverse order for pop_back)

    public:
        explicit COO(const char* value_path, const char* index_path, char mode) :
            value_file_(value_path), index_file_(index_path), mode_(mode), index_(0) {

            if (mode_ == 'w') {
                fstream_.open(value_file_, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!fstream_) throw std::runtime_error("Cannot open for write: " + value_file_);

                index_stream_.open(index_file_, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!index_stream_) throw std::runtime_error("Cannot open index for write: " + index_file_);

                // Write initial 0 to index (start of column 0)
                uint64_t zero = 0;
                index_stream_.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
            } else if (mode_ == 'r') {
                fstream_.open(value_file_, std::ios::in | std::ios::binary);
                if (!fstream_) throw std::runtime_error("File not found (read mode): " + value_file_);
            } else {
                throw std::runtime_error("mode must be 'r' or 'w'");
            }
        }

        void push(uint32_t value) {
            if (mode_ != 'w')
                throw std::runtime_error("push() not allowed in read mode");
            fstream_.write(reinterpret_cast<const char*>(&value), sizeof(value));
            index_++;
        }

        void flush_column() {
            if (mode_ != 'w')
                throw std::runtime_error("flush_column() not allowed in read mode");
            index_stream_.write(reinterpret_cast<const char*>(&index_), sizeof(index_));
        }

        void close() {
            if (fstream_.is_open()) {
                if (mode_ == 'w') {
                    fstream_.flush();
                }
                fstream_.close();
            }
            if (index_stream_.is_open()) {
                index_stream_.flush();
                index_stream_.close();
            }
        }

        void open_append() {
            if (fstream_.is_open()) {
                fstream_.close();
            }
            if (index_stream_.is_open()) {
                index_stream_.close();
            }

            // Read the last index value from the file
            std::ifstream idx_in(index_file_, std::ios::binary | std::ios::ate);
            std::streamoff file_size = idx_in.tellg();
            if (file_size >= static_cast<std::streamoff>(sizeof(uint64_t))) {
                idx_in.seekg(-static_cast<std::streamoff>(sizeof(uint64_t)), std::ios::end);
                idx_in.read(reinterpret_cast<char*>(&index_), sizeof(uint64_t));
            } else {
                index_ = 0;
            }
            idx_in.close();

            fstream_.open(value_file_, std::ios::out | std::ios::binary | std::ios::app);
            if (!fstream_) throw std::runtime_error("Cannot reopen for append: " + value_file_);

            index_stream_.open(index_file_, std::ios::out | std::ios::binary | std::ios::app);
            if (!index_stream_) throw std::runtime_error("Cannot reopen index for append: " + index_file_);

            mode_ = 'w';
        }

        void load_column(uint32_t j) const {
            if (mode_ != 'r')
                throw std::runtime_error("load_column() only valid in read mode");
            std::ifstream index_in(index_file_, std::ios::binary);
            index_in.seekg(j * sizeof(uint64_t), std::ios::beg);
            uint64_t start, end;
            index_in.read(reinterpret_cast<char*>(&start), sizeof(uint64_t));
            index_in.read(reinterpret_cast<char*>(&end), sizeof(uint64_t));

            if (start == end) {
                column_cache_.clear();
                return;
            }

            column_cache_.resize(end - start);
            fstream_.clear();
            fstream_.seekg(start * sizeof(uint32_t), std::ios::beg);
            fstream_.read(reinterpret_cast<char*>(column_cache_.data()), column_cache_.size() * sizeof(uint32_t));
            std::reverse(column_cache_.begin(), column_cache_.end());
        }

        uint32_t pop() const {
            if (mode_ != 'r')
                throw std::runtime_error("pop() not allowed in write mode");
            uint32_t value = column_cache_.back();
            column_cache_.pop_back();
            return value;
        }

        uint64_t get_index() const {
            return index_;
        }

};

} // namespace ldzip

