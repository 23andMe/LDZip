#pragma once
#include "metadata.hpp"
#include <string>
#include <vector>
#include <cstdint>

namespace ldzip {

void concat_ldzip(const std::vector<std::string> &prefixes,
                        const std::string &out_prefix);

// Concatenate overlapping chunks (tabularPlink / FULL format only).
// Overlap between consecutive chunks is auto-detected from the vars files.
// For each column in the overlapping region the above-diagonal entries
// (row <= column) are taken from the earlier chunk and the below-diagonal
// entries (row > column) are taken from the later chunk.
void concat_ldzip_overlapping(const std::vector<std::string> &prefixes,
                              const std::string &out_prefix);

} // namespace ldzip
