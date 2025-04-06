#ifndef __UTILITIES__
#define __UTILITIES__

#include <Rcpp.h>
#include <algorithm>  // for std::copy
#include <HydroGallery.h>
using namespace Rcpp;
using namespace HydroGallery;

void save_wgmat(const NumericMatrix& matrix, const std::string& filename);
NumericMatrix load_wgmat(const std::string& filename);
void bind_wgmat(const StringVector& input_files, const std::string& output_file);

#endif // __UTILITIES__
