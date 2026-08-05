#ifndef FIND_TAG_VARIANTS_HPP
#define FIND_TAG_VARIANTS_HPP

#include <string>
#include <vector>

namespace ldzip {
    void find_tag_variants(
        const std::string& input_prefix,
        const std::string& variant_file,
        const std::string& output_file,
        double threshold,
        const std::string& stat
    );
}

#endif // FIND_TAG_VARIANTS_HPP
