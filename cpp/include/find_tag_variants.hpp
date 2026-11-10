#pragma once
#include <string>

namespace ldzip {

void find_tag_variants(const std::string &input_prefix,
                       const std::string &variant_file,
                       const std::string &output_file,
                       double threshold,
                       const std::string &stat);

} // namespace ldzip
