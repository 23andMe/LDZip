#include "CLI11.hpp"
#include <filesystem>
#include "binary_compressor.hpp"
#include "tabular_compressor.hpp"
#include "binary_decompressor.hpp"
#include "concat.hpp"
#include "filter.hpp"


int main(int argc, char** argv) {
    CLI::App app{"ldzip - compress/decompress LD matrices"};

    // Subcommand: compress
    const std::vector<std::string> column_names = {"PHASED_R", "UNPHASED_R", "PHASED_R2", "UNPHASED_R2", "D", "DPRIME"};
    std::string ld_file, snp_file, output_prefix, keep_file, format = "full", type="binary", min_col, range_str="";
    std::vector<std::string> columns{"PHASED_R"};  // from CLI
    std::vector<ldzip::Stat> stats;
    ldzip::Stat min_stat;
    int bits = 8;
    float min = 0.01;
    static const std::map<std::string, ldzip::Stat> col_map {
        {"PHASED_R",    ldzip::Stat::PHASED_R},
        {"UNPHASED_R",  ldzip::Stat::UNPHASED_R},
        {"PHASED_R2",   ldzip::Stat::PHASED_R2},
        {"UNPHASED_R2", ldzip::Stat::UNPHASED_R2},
        {"D",           ldzip::Stat::D},
        {"DPRIME",      ldzip::Stat::DPRIME},
    };

    auto compress = app.add_subcommand("compress", "Compress a binary LD matrix");
    auto filter = app.add_subcommand("filter", "Filter a .ldzip binary matrix");
    auto decompress = app.add_subcommand("decompress", "Decompress a .ldzip binary matrix");
    auto plink_bin = compress->add_subcommand("plinkSquare","Plink Binary Input");
    auto plink_tabular = compress->add_subcommand("plinkTabular","Plink Tabular Input");

    plink_bin->add_option("-l,--ld_file", ld_file, "Input Plink Binary Matrix")->required()->check(CLI::ExistingFile);
    plink_bin->add_option("-s,--snp_file", snp_file, "Variant list file")->required()->check(CLI::ExistingFile);
    plink_bin->add_option("-o,--output_prefix", output_prefix, "Output prefix")->required();
    plink_bin->add_option("-b,--bits", bits, "Number of bits to quanitze LD values (99 - float)")->check(CLI::IsMember({8, 16, 32, 99}))->default_val(std::to_string(bits));
    plink_bin->add_option("-m,--min", min, "Minimum absolute threshold")->check(CLI::NonNegativeNumber)->default_val(std::to_string(min));
    plink_bin->add_option("-t,--type", columns, "Type of LD available in binary matrix. ")->check(CLI::IsMember(column_names, CLI::ignore_case))->default_val(columns[0])->expected(1); 
 
    
    plink_tabular->add_option("-l,--ld_file", ld_file, "Input LD file")->required()->check(CLI::ExistingFile);
    plink_tabular->add_option("-s,--snp_file", snp_file, "SNP list file")->required()->check(CLI::ExistingFile);
    plink_tabular->add_option("-o,--output_prefix", output_prefix, "Output prefix")->required();
    plink_tabular->add_option("-b,--bits", bits, "Bits (8,16,32,99)")->check(CLI::IsMember({8, 16, 32, 99}))->default_val(std::to_string(bits));
    plink_tabular->add_option("-c,--columns", columns, "LD columns to use from plink tabular file (comma separated)")->check(CLI::IsMember(column_names, CLI::ignore_case))->default_val(columns[0])->delimiter(','); 
    plink_tabular->add_option("-k,--min_col", min_col, "LD column to apply minimum absolute thresholding")->check(CLI::IsMember(column_names, CLI::ignore_case))->default_val(columns[0])->expected(1); 
    plink_tabular->add_option("-m,--min", min, "Minimum absolute threshold")->check(CLI::NonNegativeNumber)->default_val(std::to_string(min));

    // Subcommand: filter
    std::string input_prefix;
    filter->add_option("-i,--input_prefix",     input_prefix,   "Input .ldzip prefix (expects .x.bin, .i.bin, .io.bin, .io.index, .p.bin, .meta.json)")->required();
    filter->add_option("-o,--output_prefix",    output_prefix,  "Output prefix for filtered matrix")->required();
    auto range_opt = filter->add_option("-r,--range", range_str, "Range of indices (0-based START-END)");
    auto keep_opt  = filter->add_option("-k,--keep",  keep_file, "File with list of indices (0-based)");
    range_opt->excludes(keep_opt);

    // Subcommand: decompress
    decompress->add_option("-i,--input_prefix",     input_prefix,   "Input .ldzip prefix (expects .x.bin, .i.bin, .io.bin, .io.index, .p.bin, .meta.json)")->required();
    decompress->add_option("-t,--type",             type,           "Output Type (tabular or binary)")->default_val(type);
    decompress->add_option("-o,--output_prefix",    output_prefix,  "Output path for decompressed matrix")->required();

    // Subcommand: concat
    std::vector<std::string> input_chunks;
    std::string concat_output;
    bool overlapping = false;
    auto concat = app.add_subcommand("concat", "Concatenate multiple .ldzip chunks");
    concat->add_option("-i,--inputs", input_chunks, "List of input chunk prefixes")->required()->expected(-1);
    concat->add_option("-o,--output_prefix", concat_output, "Output prefix for concatenated .ldzip")->required();
    concat->add_flag("--overlapping", overlapping, "Input chunks have overlapping variant regions");

    CLI11_PARSE(app, argc, argv);

    stats.reserve(columns.size());
    for (auto &s : columns) {
        stats.push_back(ldzip::parse_stat(s));
    }
    min_stat = ldzip::parse_stat(min_col);
    if (app.get_subcommands().empty()) {
        std::cout << app.help() << std::endl;
        std::exit(1);
    }

    if (compress->parsed()) {

        if (compress->get_subcommands().empty()) {
            std::cout << compress->help() << std::endl;
            std::exit(1);
        }

        if(plink_bin->parsed()){
            
            ldzip::compress_binary_matrix(ld_file, snp_file, output_prefix, bits, min, format, stats[0]);

        }else if(plink_tabular->parsed()){

            ldzip::compress_tabular_file(ld_file, snp_file, output_prefix, bits, min, format, min_stat);
            
        }
    } else if(filter->parsed()) {
        
        std::vector<size_t> indices;
        if (!range_str.empty()) {
            int start, end;
            if (std::sscanf(range_str.c_str(), "%d-%d", &start, &end) != 2 || start > end)
                throw std::runtime_error("Invalid range format: expected START-END (0-based)");

            indices.resize(end - start + 1);
            for (int i = 0; i <= end - start; ++i)
                indices[i] = start + i;

        } else if (!keep_file.empty()) {
            std::ifstream in(keep_file);
            if (!in)
                throw std::runtime_error("Failed to open keep file: " + keep_file);

            std::string line;
            int line_no = 0, idx;
            char extra;

            while (std::getline(in, line)) {
                ++line_no;
                std::istringstream iss(line);
                if (!(iss >> idx) || (iss >> extra)) 
                    throw std::runtime_error("Invalid line " + std::to_string(line_no) + " in keep file (must contain exactly one integer)");
                if (idx < 0)
                    throw std::runtime_error("Negative index at line " + std::to_string(line_no));
                indices.push_back(idx);
            }
        }

        ldzip::filter_ldzip(input_prefix, output_prefix, indices);

    } else if (decompress->parsed()) {

        ldzip::decompress_ldzip(input_prefix, output_prefix, type);

    } else if (concat->parsed()) {

        ldzip::concat_ldzip(input_chunks, concat_output, overlapping);
    }


    return 0;
}
