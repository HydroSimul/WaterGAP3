#ifndef __UTILITIES__
#define __UTILITIES__

#include <Rcpp.h>
#include <algorithm>  // for std::copy
using namespace Rcpp;

NumericVector vecpow(NumericVector base, NumericVector exp);
NumericVector vecpow10(NumericVector exp);
double sum_product(NumericVector lhs, NumericVector rhs);
void resetVector(Rcpp::NumericVector& x);
NumericVector subset_get(NumericVector num_Data, IntegerVector int_Index);
LogicalVector subset_get_logical(LogicalVector data, IntegerVector int_Index);
void subset_put(NumericVector &num_Data, IntegerVector int_Index, NumericVector num_DataPut);
IntegerVector c_int(IntegerVector vec_A, IntegerVector vec_B);
NumericVector c_num(NumericVector vec_A, NumericVector vec_B);
#endif // __UTILITIES__