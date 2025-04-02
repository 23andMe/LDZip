#include <Rcpp.h>
#include "ldzipmatrix.hpp"

using namespace Rcpp;


// [[Rcpp::export]]
SEXP LDZipMatrix_rcpp(std::string prefix) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(new ldzip::LDZipMatrix(prefix), true);
    return ptr;
}

// [[Rcpp::export]]
size_t LDZipMatrix_nrows_rcpp(SEXP ld) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    return ptr->nrows();
}

// [[Rcpp::export]]
size_t LDZipMatrix_ncols_rcpp(SEXP ld) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    return ptr->ncols();
}

// [[Rcpp::export]]
uint64_t LDZipMatrix_nnz_rcpp(SEXP ld) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    return ptr->nnz();
}

// [[Rcpp::export]]
int LDZipMatrix_bits_rcpp(SEXP ld) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    return ptr->bits();  // int version
}

// [[Rcpp::export]]
Rcpp::CharacterVector LDZipMatrix_allStats_rcpp() {
    auto arr = ldzip::All_Stats();
    Rcpp::CharacterVector out(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        out[i] = ldzip::stat_to_string(arr[i]); // assuming you have a helper
    }
    return out;
}

// [[Rcpp::export]]
std::string LDZipMatrix_format_rcpp(SEXP ld) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    return (ptr->formatStr());
}

// [[Rcpp::export]]
double LDZipMatrix_getValue_rcpp(SEXP ld, int row, int col, std::string type) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    return ptr->getValue(row, col, type);
}


// [[Rcpp::export]]
bool LDZipMatrix_hasStatString_rcpp(SEXP ld, std::string type) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    return ptr->has_stat(type);
}

// [[Rcpp::export]]
Rcpp::NumericVector LDZipMatrix_get_column_range_rcpp(
                                        SEXP ld,
                                        size_t column,
                                        size_t row_start,
                                        size_t row_end, 
                                        std::string type) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    Rcpp::NumericVector out(row_end - row_start + 1);
    ptr->getColumn(column, row_start, row_end, REAL(out), type);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector LDZipMatrix_get_column_indices_rcpp(
                                        SEXP ld,
                                        size_t column,
                                        Rcpp::IntegerVector rows, 
                                        std::string type) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    Rcpp::NumericVector out(rows.size());
    ptr->getColumn(column, rows.begin(), rows.size(), REAL(out), type);
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerVector LDZipMatrix_get_i_rcpp(
    
                                        SEXP ld, 
                                        size_t column) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    std::vector<uint32_t> idx = ptr->get_i(static_cast<uint32_t>(column));
    return Rcpp::wrap(idx);
}


// [[Rcpp::export]]
Rcpp::IntegerVector LDZipMatrix_get_neighbors_rcpp(
                                        SEXP ld, 
                                        size_t column,
                                        double abs_threshold,
                                        std::string type) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    std::vector<uint32_t> idx;
    idx = ptr->get_neighbors(static_cast<uint32_t>(column), (abs_threshold), type);
    return Rcpp::wrap(idx);
}

// [[Rcpp::export]]
std::string LDZipMatrix_get_prefix_rcpp(SEXP ld) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    return ptr->prefix();
}


// [[Rcpp::export]]
Rcpp::NumericMatrix LDZipMatrix_get_submatrix_range_rcpp(
                                        SEXP ld,
                                        size_t col_start,
                                        size_t col_end,
                                        size_t row_start,
                                        size_t row_end,
                                        std::string type) {
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    Rcpp::NumericMatrix out(row_end - row_start + 1, col_end - col_start + 1); 
    ptr->getSubMatrix(col_start, col_end, row_start, row_end, REAL(out), type);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LDZipMatrix_get_submatrix_indices_rcpp(
    SEXP ld,
    Rcpp::IntegerVector cols,
    Rcpp::IntegerVector rows,
    std::string type) 
{
    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);  
    Rcpp::NumericMatrix out(rows.size(), cols.size());

    // pass raw int* directly (no copies)
    ptr->getSubMatrix(cols.begin(), cols.size(),
                      rows.begin(), rows.size(),
                      out.begin(), type);

    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector LDZipMatrix_get_pairswise_rcpp(
    SEXP ld,
    Rcpp::IntegerVector cols,
    Rcpp::IntegerVector rows,
    std::string type) {

    Rcpp::XPtr<ldzip::LDZipMatrix> ptr(ld);
    Rcpp::NumericVector out(rows.size());
    ptr->getPairwise(cols.begin(), rows.begin(), rows.size(), REAL(out), type);
    return out;
}