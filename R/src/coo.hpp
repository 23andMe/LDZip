#include <fstream>
#include <cstdint>
#include <unordered_map>
#include <optional>
#include <stdexcept>
#include <string>
#include <iostream>

namespace ldzip
{
class COO {
    public:
        struct Entry {
            uint32_t i;
            uint32_t j;
            uint32_t value;
        };

        using Key = uint64_t;

    private:
        std::fstream fstream_;
        std::unordered_map<Key, uint32_t> hash_;
        char mode_;

    public:
        COO() = default;

        explicit COO(const char* path, char mode) : mode_(mode) {
            if (mode_ == 'w') {
                fstream_.open(path, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!fstream_) throw std::runtime_error("Cannot open for write: " + std::string(path));
            } else if (mode_ == 'r') {
                fstream_.open(path, std::ios::in | std::ios::binary);
                if (!fstream_) throw std::runtime_error("File not found (read mode): " + std::string(path));
            } else {
                throw std::runtime_error("mode must be 'r' or 'w'");
            }
        }   

        void push(uint32_t i, uint32_t j, uint32_t value) {
            if (mode_ != 'w')
                throw std::runtime_error("push() not allowed in read mode");
            Entry e{i, j, value};
            fstream_.seekp(0, std::ios::end);
            fstream_.write(reinterpret_cast<const char*>(&e), sizeof(e));
            fstream_.flush();
        }

        static Key pack(uint32_t i, uint32_t j) {
            return (uint64_t(i) << 32) | j;
        }

        void load() {
            hash_.clear();
            fstream_.clear();
            fstream_.seekg(0, std::ios::end);

            std::streamoff bytes = fstream_.tellg();
            if (bytes < 0 || bytes % sizeof(Entry) != 0)
                throw std::runtime_error("Invalid COO file size");
            size_t nnz = bytes / sizeof(Entry);
            hash_.reserve(nnz);

            fstream_.seekg(0, std::ios::beg);
            for (size_t k = 0; k < nnz; ++k) {
                Entry e;
                fstream_.read(reinterpret_cast<char*>(&e), sizeof(e));
                if (!fstream_)
                    throw std::runtime_error("Error reading COO entry");
                hash_.emplace(pack(e.i, e.j), e.value);
            }
        }

        uint32_t get(uint32_t i, uint32_t j) const {
            if (mode_ != 'r')
                throw std::runtime_error("get() not allowed in write mode");
            auto it = hash_.find(pack(i, j));
            if (it == hash_.end())
                throw std::out_of_range("COO entry not found");
            return it->second;
        }

};

} // namespace ldzip

