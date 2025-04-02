#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include "metadata.hpp"

namespace snp_util {

int count_lines(const std::string& filename);

int countFields(const std::string& filename);

std::vector<std::string> getFields(const std::string& filename, bool no_header = false);
void write_snp_subset(const std::string& inFile,
                      const std::string& outFile,
                      const std::vector<size_t>& keep = {});

void get_snp_hash   (const std::string& filename, 
                    std::unordered_map<std::string, uint32_t> &rsid_to_index);

std::pair<std::vector<ldzip::Stat>, std::unordered_map<int, ldzip::Stat>>
    get_stats_from_fields(const std::vector<std::string>& fields);
    

}
