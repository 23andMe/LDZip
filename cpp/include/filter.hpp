#pragma once
#include "metadata.hpp"
#include <string>
#include <vector>

namespace ldzip {

void filter_ldzip(const std::string &in_prefix,
                const std::string &out_prefix,
                std::vector<size_t> indices);

} // namespace ldzip
