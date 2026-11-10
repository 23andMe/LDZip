#pragma once
#include <string>

namespace ldzip {

void ld_score(const std::string &input_prefix,
              const std::string &output_file,
              size_t window_kb,
              float threshold,
              const std::string &stat);

} // namespace ldzip
