#ifndef LD_PRUNING_HPP
#define LD_PRUNING_HPP

#include <string>
#include <cstdint>

namespace ldzip {

    void ld_pruning(
        const std::string& input_prefix,
        const std::string& output_prefix,
        double threshold,
        size_t window_kb,
        const std::string& stat
    );

}

#endif // LD_PRUNING_HPP
