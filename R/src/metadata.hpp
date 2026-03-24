#pragma once
#include "json.hpp"
#include <string>
#include <stdexcept>
#include <array>
#include <fstream>

// A wrapper around std::array that lets us use an enum as the index
template <typename T, typename Enum, size_t N = static_cast<size_t>(Enum::NTYPES)>
struct EnumArray {
    std::array<T, N> data;

    // Default constructor (all zeros)
    EnumArray() : data{} {}

    // Fill constructor
    explicit EnumArray(const T& value) {
        data.fill(value);
    }
    
    // Overload operator[] so we can use Enum directly
    T& operator[](Enum e) { 
        return data[static_cast<size_t>(e)]; 
    }
    const T& operator[](Enum e) const { 
        return data[static_cast<size_t>(e)]; 
    }

    // Allow iteration in range-based for loops
    auto begin() { return data.begin(); }
    auto end()   { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end()   const { return data.end(); }
};

namespace ldzip {

    // Defines storage format for LD matrix
    enum class MatrixFormat { FULL, UPPER };
    enum class Bits : uint16_t {
        B8  = 8,
        B16 = 16,
        B32 = 32,
        B99 = 99
    };
    enum class Stat : uint8_t { 
        PHASED_R = 0, 
        UNPHASED_R = 1, 
        PHASED_R2 = 2,
        UNPHASED_R2 = 3,
        D = 4, 
        DPRIME = 5, 
        COUNT = 6,
        NTYPES = 7 };
    constexpr auto All_Stats() {
        std::array<Stat, static_cast<size_t>(Stat::NTYPES)> arr{};
        for (size_t i = 0; i < arr.size(); ++i) {
            arr[i] = static_cast<Stat>(i);
        }
        return arr;
    }

    // Convert MatrixFormat enum to string ("full", "upper")
    inline std::string format_to_string(MatrixFormat fmt) {
        switch (fmt) {
            case MatrixFormat::FULL: return "full";
            case MatrixFormat::UPPER: return "upper";
            default: throw std::invalid_argument("Unknown MatrixFormat");
        }
    }

    // Convert MatrixFormat string to enum
    inline MatrixFormat parse_format(const std::string& str) {
        if (str == "full") return MatrixFormat::FULL;
        if (str == "upper") return MatrixFormat::UPPER;
        throw std::invalid_argument("Unsupported Matrix format: " + str);
    }

    // Convert Bits enum to integer (8, 16, 32, 99)
    inline uint16_t bits_to_int(Bits b) {
        return static_cast<uint16_t>(b);
    }

    // Convert Bits integer to enum
    inline Bits parse_bits(uint16_t value) {
        switch (value) {
            case 8:  return Bits::B8;
            case 16: return Bits::B16;
            case 32: return Bits::B32;
            case 99: return Bits::B99;
        }
        throw std::invalid_argument("Unsupported Bits value: " + std::to_string(value));
    }

    // Convert Stat enum to string ("r", "r2", "d", "dprime", "n")
    inline std::string stat_to_string(Stat s) {
        switch (s) {
            case Stat::PHASED_R: return "PHASED_R";
            case Stat::UNPHASED_R: return "UNPHASED_R";
            case Stat::PHASED_R2: return "PHASED_R2";
            case Stat::UNPHASED_R2: return "UNPHASED_R2";
            case Stat::D: return "D";
            case Stat::DPRIME: return "DPRIME";
            case Stat::COUNT: return "N";
            default: throw std::invalid_argument("Unknown Stat");
        }
    }

    // Convert Stat string to enum
    inline Stat parse_stat(const std::string& str) {
        if (str == "PHASED_R") return Stat::PHASED_R;
        if (str == "PHASED_R2") return Stat::PHASED_R2;
        if (str == "UNPHASED_R") return Stat::UNPHASED_R;
        if (str == "UNPHASED_R2") return Stat::UNPHASED_R2;
        if (str == "D") return Stat::D;
        if (str == "DPRIME") return Stat::DPRIME;
        if (str == "N") return Stat::COUNT;
        throw std::invalid_argument("Unsupported Stat/Type format: " + str);
    }

    class Variant {
        public:
            Variant() = default;

            Variant(std::string id, std::string ref, std::string alt)
                : id_(std::move(id)), ref_(std::move(ref)), alt_(std::move(alt)) {}

            // --- Getters ---
            const std::string& id()  const { return id_; }
            const std::string& ref() const { return ref_; }
            const std::string& alt() const { return alt_; }

            // --- Setters ---
            void set_id(const std::string& id)   { id_ = id; }
            void set_ref(const std::string& ref) { ref_ = ref; }
            void set_alt(const std::string& alt) { alt_ = alt; }

        private:
            std::string id_;
            std::string ref_;
            std::string alt_;
    };

    // Metadata for the LDZip matrix
    struct MetaInfo {
        size_t rows = 0;
        size_t cols = 0;
        uint64_t nnz = 0;
        Bits bits = Bits::B8;
        MatrixFormat format = MatrixFormat::UPPER;
        EnumArray<bool, Stat> has_stat{};
        std::string version = "";

        MetaInfo() = default;

        MetaInfo(size_t r, size_t c, uint64_t n, int bit_val, MatrixFormat fmt, std::string ver)
            : rows(r), cols(c), nnz(n), bits(parse_bits(bit_val)), format(fmt), version(ver) {}
    };

    inline void write_metadata_json(const std::string& path, const MetaInfo& meta) {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot open output metadata file : " + path);

        nlohmann::json j;
        j["dim"]["rows"] = meta.rows;
        j["dim"]["cols"] = meta.cols;
        j["nnz"] = meta.nnz;
        j["bits"] = bits_to_int(meta.bits);
        j["format"] = format_to_string(meta.format);
        j["version"] = meta.version;

        nlohmann::json stats_json = nlohmann::json::array();
        for (int i = 0; i < static_cast<int>(Stat::NTYPES); ++i) {
            Stat s = static_cast<Stat>(i);
            if (meta.has_stat[s]) {
                stats_json.push_back(stat_to_string(s));
            }
        }
        j["stats"] = stats_json;

        out << j.dump(2) << std::endl;  // pretty print with 2-space indent
    }
                   
    inline MetaInfo read_metadata_json(const std::string& path) {
        std::ifstream in(path);
        if (!in.is_open()) throw std::runtime_error("Could not open input metadata file: " + path);

        nlohmann::json j;
        in >> j;

        MetaInfo meta;
        meta.rows = j["dim"]["rows"];
        meta.cols = j["dim"]["cols"];
        meta.bits = parse_bits(j["bits"]);
        meta.nnz = j["nnz"];
        meta.format = parse_format(j["format"]);
        meta.version = "1.1";
        if (j.contains("version")) {
            meta.version = j["version"].get<std::string>().c_str();
        }
        for (const auto& stat_str : j["stats"]) {
            Stat s = parse_stat(stat_str.get<std::string>());
            meta.has_stat[s] = true;
        }
        return meta;
    }

}
