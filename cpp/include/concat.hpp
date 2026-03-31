#pragma once
#include "metadata.hpp"
#include <string>
#include <vector>
#include <cstdint>

namespace ldzip {

void concat_ldzip(const std::vector<std::string> &prefixes,
                  const std::string &out_prefix,
                  bool overlapping = false);

} // namespace ldzip
