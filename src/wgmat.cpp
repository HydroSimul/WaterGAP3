#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <vector>

// [[Rcpp::interfaces(r, cpp)]]

//' Save a NumericMatrix to a binary file
//'
//' This function saves a NumericMatrix to a binary file, with the matrix
//' data stored in a row-major order (transposed).
//' @rdname wgmat
//' @param matrix The NumericMatrix to be saved.
//' @param filename The name of the binary file to save the matrix to.
//'
//' @export
// [[Rcpp::export]]
void save_wgmat(const NumericMatrix& matrix, const std::string& filename) {
 int rows = matrix.nrow();
 int cols = matrix.ncol();
 std::ofstream file(filename, std::ios::binary);
 if (!file.is_open()) {
   stop("Failed to open file for writing: " + filename);
 }

 // Write dimensions
 file.write(reinterpret_cast<const char*>(&rows), sizeof(int));
 file.write(reinterpret_cast<const char*>(&cols), sizeof(int));

 // Write matrix data in row-major order (transposed)
 std::vector<float> buffer(rows * cols);
 for (int i = 0; i < rows; ++i) {
   for (int j = 0; j < cols; ++j) {
     buffer[j + i * cols] = static_cast<float>(matrix(i, j));  // Transpose during saving
   }
 }

 file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(float));
 file.close();
}

//' Load a NumericMatrix from a binary file
//'
//' This function loads a NumericMatrix from a binary file that was saved
//' using the `save_wgmat` function. The matrix is reconstructed in R's
//' column-major order.
//' @name wgmat
//' @param filename The name of the binary file to load the matrix from.
//'
//' @return A NumericMatrix containing the data from the binary file.
//' @export
// [[Rcpp::export]]
NumericMatrix load_wgmat(const std::string& filename) {
 int rows, cols;
 std::ifstream file(filename, std::ios::binary);
 if (!file.is_open()) {
   stop("Failed to open file for reading: " + filename);
 }

 // Read dimensions
 file.read(reinterpret_cast<char*>(&rows), sizeof(int));
 file.read(reinterpret_cast<char*>(&cols), sizeof(int));

 // Read matrix data into buffer
 std::vector<float> buffer(rows * cols);
 file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(float));
 file.close();

 // Reconstruct matrix in R's column-major order
 NumericMatrix matrix(rows, cols);
 for (int i = 0; i < rows; ++i) {
   for (int j = 0; j < cols; ++j) {
     matrix(i, j) = static_cast<double>(buffer[j + i * cols]);  // Re-transpose after loading
   }
 }

 return matrix;
}


//' Bind multiple matrices into a single binary file
//'
//' This function binds several matrices from binary files into a single
//' matrix and saves it to a new binary file. All input matrices must have
//' the same number of columns. The function validates the dimensions of
//' the matrices before combining them.
//' @name wgmat
//' @param input_files A vector of file paths to the input binary files.
//' @param output_file The file path to save the combined matrix.
//'
//' @export
// [[Rcpp::export]]
void bind_wgmat(const StringVector& input_files, const std::string& output_file) {

 // Check if we have at least one input file
 if (input_files.size() == 0) {
   stop("No input files provided");
 }

 int total_rows = 0;
 int total_cols = 0;
 int expected_cols = 0;

 // First pass: load all files and validate dimensions
 std::vector<std::vector<float>> all_buffers;

 for (int f = 0; f < input_files.size(); f++) {
   std::string filename = as<std::string>(input_files[f]);

   // Open the file to read the matrix
   int rows, cols;
   std::ifstream file(filename, std::ios::binary);
   if (!file.is_open()) {
     stop("Failed to open file for reading: " + filename);
   }

   // Read the dimensions (rows and columns)
   file.read(reinterpret_cast<char*>(&rows), sizeof(int));
   file.read(reinterpret_cast<char*>(&cols), sizeof(int));

   // Validate that the columns match the expected number of columns
   if (f == 0) {
     expected_cols = cols;
     total_cols = expected_cols;
   } else {
     if (cols != expected_cols) {
       stop("Column mismatch in file " + filename + ": expected " +
         std::to_string(expected_cols) + " columns, found " +
         std::to_string(cols));
     }
   }

   // Read the matrix data into a buffer
   std::vector<float> buffer(rows * cols);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(float));
   file.close();

   // Store the buffer for this matrix
   total_rows += rows;
   all_buffers.push_back(buffer);
 }

 // Create output file
 std::ofstream outfile(output_file, std::ios::binary);
 if (!outfile.is_open()) {
   stop("Failed to open output file for writing: " + output_file);
 }

 // Write the header with final dimensions (total_rows, total_cols)
 outfile.write(reinterpret_cast<const char*>(&total_rows), sizeof(int));
 outfile.write(reinterpret_cast<const char*>(&total_cols), sizeof(int));

 // Write the combined matrix data row by row
 for (const auto& buffer : all_buffers) {
   outfile.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(float));
 }

 // Close the output file
 outfile.close();
}
