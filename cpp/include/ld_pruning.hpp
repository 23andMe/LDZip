#pragma once
#include <string>

namespace ldzip {

void ld_pruning(const std::string &input_prefix,
                const std::string &output_prefix,
                double threshold,
                size_t window_kb,
                const std::string &stat);

} // namespace ldzip
