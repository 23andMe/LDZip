#include "snp_util.hpp"
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>

namespace snp_util {

int count_lines(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        ++count;
    }
    return count;
}

int countFields(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);
    std::string line;
    while (std::getline(in, line)) 
        if (!line.empty() && line[0] != '#') {
            return std::count(line.begin(), line.end(), '\t') + 1;
        }
    throw std::runtime_error("No single-hash header line found in file: " + filename);
}

std::vector<std::string> getFields(const std::string& filename, bool no_header) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && (no_header || (line[0] == '#' && (line.size() == 1 || line[1] != '#')))) {
            if (!no_header) line.erase(0, 1);
            std::istringstream iss(line);
            std::vector<std::string> fields;
            std::string field;
            while (std::getline(iss, field, '\t')) fields.push_back(field);
            return fields;
        }
    }
    throw std::runtime_error("No header line found in file: " + filename);
}



void write_snp_subset(const std::string& inFile,
                      const std::string& outFile,
                      const std::vector<size_t>& keep) {
    std::ifstream in(inFile);
    if (!in) throw std::runtime_error("Cannot open file: " + inFile);
    std::ofstream out(outFile);
    if (!out) throw std::runtime_error("Cannot open file: " + outFile);

    std::string line;
    bool has_header = false;

    // detect header
    std::streampos start = in.tellg();
    if (std::getline(in, line)) {
        auto it = std::find_if_not(line.begin(), line.end(), ::isspace);
        has_header = (it != line.end() && !std::isdigit(*it));
        if (has_header) out << line << '\n';
        else in.seekg(start);
    }

    // if keep is empty → copy everything
    if (keep.empty()) {
        while (std::getline(in, line)) 
            out << line << '\n';
        return;
    }

    size_t data_idx = 0;
    size_t keep_idx = 0;

    while (std::getline(in, line) && keep_idx < keep.size()) {
        if (data_idx == keep[keep_idx]) {
            out << line << '\n';
            ++keep_idx;
        }
        ++data_idx;
    }
}

void get_snp_hash(const std::string& filename, std::unordered_map<std::string, uint32_t> &rsid_to_index) {

    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open SNP file: " + filename);

    size_t NoFields = 5;
    std::string line;
    uint32_t index = 0;
    std::vector<std::string> fields = getFields(filename);
        if (fields.size() < NoFields) throw std::runtime_error("Expected at least 5 fields, found " + std::to_string(fields.size()) + " in " + filename);
    const std::vector<std::string> required = {"CHROM","POS","ID","REF","ALT"};
    for (size_t i = 0; i < required.size(); ++i)
        if (fields[i] != required[i]) throw std::runtime_error("Expected " + required[i] + " at col " + std::to_string(i+1) + ", found " + fields[i] + " in " + filename);

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        if ((index+1) % 1000000 == 0) std::cout << " Hashed " << index+1 << " entries" << std::endl;
        for (size_t i = 0; i < NoFields; ++i) {
            if (!std::getline(ss, fields[i], '\t'))
                throw std::runtime_error("Line has fewer than 5 columns: " + line);
        }

        const std::string& rsid = fields[2];
        if (rsid_to_index.find(rsid) != rsid_to_index.end()) {
            throw std::runtime_error(
                "Duplicate variant ID '" + rsid + "' found in " + filename +
                " : first at non-header line " + std::to_string(rsid_to_index[rsid] + 1) +
                ", again at non-header line " + std::to_string(index + 1)
            );
        }

        rsid_to_index[rsid] = index++;
    }

}


std::pair<std::vector<ldzip::Stat>, std::unordered_map<int, ldzip::Stat>>
    get_stats_from_fields(const std::vector<std::string>& fields) {
    
    std::vector<ldzip::Stat> stats_available;
    std::unordered_map<int, ldzip::Stat> col_to_stat;

    for (size_t col = 0; col < fields.size(); ++col) {
        try {
            ldzip::Stat stat = ldzip::parse_stat(fields[col]);
            stats_available.push_back(stat);
            col_to_stat[static_cast<int>(col)] = stat;
        } catch (const std::invalid_argument&) {
            // ignore non-stat columns
        }
    }

    return {stats_available, col_to_stat};
}



} // snp_util
