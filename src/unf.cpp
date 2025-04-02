#include <Rcpp.h>
#include <string>
#include <fstream>
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]


List analyse_file_unf(std::string fn_UNF) {
  // Extract the last character and convert to integer
  int mark_UNF = fn_UNF[fn_UNF.length() - 1] - '0';

  // Create vector of type names
  std::vector<std::string> type_vec = {"numeric", "int", "int", "", "integer", "", "", "", "double"};
  std::string type_UNF = (mark_UNF >= 0 && mark_UNF < (int)type_vec.size()) ? type_vec[mark_UNF] : "";

  // Calculate number of bytes
  int n_Byte = mark_UNF > 0 ? mark_UNF : 4;

  return List::create(
    Named("type_UNF") = type_UNF,
    Named("n_Byte") = n_Byte
  );
}

//' read_UNF - Read data from UNF file
//' @name unf
//' @param fn_UNF A string, name of the UNF file
//' @return A numeric or integer vector by reading the UNF file
//' @export
// [[Rcpp::export]]
SEXP read_unf(std::string fn_UNF) {
  // Get file information
  List info_UNF = analyse_file_unf(fn_UNF);
  std::string type_UNF = info_UNF["type_UNF"];
  int n_Byte = info_UNF["n_Byte"];

  // Get file size
  std::ifstream file(fn_UNF, std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    stop("Cannot open file: " + fn_UNF);
  }

  std::streamsize n_FileSize = file.tellg();
  int n_DataSize = n_FileSize / n_Byte;
  file.seekg(0, std::ios::beg);

  // Improved endian swap functions for signed integers
  auto swap_endian_int32 = [](int32_t x) {
    return ((x & 0x000000FF) << 24) |
      ((x & 0x0000FF00) << 8) |
      ((x & 0x00FF0000) >> 8) |
      ((x & 0xFF000000) >> 24);
  };

  auto swap_endian_int16 = [](int16_t x) {
    return ((x & 0x00FF) << 8) |
      ((x & 0xFF00) >> 8);
  };

  auto swap_endian_double = [](double x) {
    double result;
    char* src = (char*)&x;
    char* dst = (char*)&result;

    for (int i = 0; i < 8; i++) {
      dst[i] = src[7-i];
    }

    return result;
  };

  auto swap_endian_float = [](float x) {
    float result;
    char* src = (char*)&x;
    char* dst = (char*)&result;

    for (int i = 0; i < 4; i++) {
      dst[i] = src[3-i];
    }

    return result;
  };

  // Read based on data type
  if (type_UNF == "numeric") {
    // For UNF0 (4-byte numeric/float)
    std::vector<float> buffer(n_DataSize);
    file.read(reinterpret_cast<char*>(&buffer[0]), n_FileSize);

    // Convert to NumericVector for R with endian swap
    NumericVector result(n_DataSize);
    for (int i = 0; i < n_DataSize; i++) {
      result[i] = swap_endian_float(buffer[i]);
    }
    return result;
  }
  else if (type_UNF == "double") {
    // For UNF8 (8-byte numeric/double)
    std::vector<double> buffer(n_DataSize);
    file.read(reinterpret_cast<char*>(&buffer[0]), n_FileSize);

    // Convert to NumericVector for R with endian swap
    NumericVector result(n_DataSize);
    for (int i = 0; i < n_DataSize; i++) {
      result[i] = swap_endian_double(buffer[i]);
    }
    return result;
  }
  else if (type_UNF == "int") {
    if (n_Byte == 1) {
      // For UNF1 (1-byte integer - no endian swap needed)
      std::vector<int8_t> buffer(n_DataSize);
      file.read(reinterpret_cast<char*>(&buffer[0]), n_FileSize);

      IntegerVector result(n_DataSize);
      for (int i = 0; i < n_DataSize; i++) {
        result[i] = buffer[i];  // No need for byte swapping
      }
      return result;
    }
    else if (n_Byte == 2) {
      // For UNF2 (2-byte integer)
      // Read as raw bytes first
      std::vector<unsigned char> raw_buffer(n_FileSize);
      file.read(reinterpret_cast<char*>(&raw_buffer[0]), n_FileSize);

      IntegerVector result(n_DataSize);
      for (int i = 0; i < n_DataSize; i++) {
        // Manually construct the int16_t with correct byte order (big endian)
        int16_t value = (raw_buffer[i*2] << 8) | raw_buffer[i*2 + 1];
        result[i] = value;
      }
      return result;
    }
  }
  else if (type_UNF == "integer") {
    // For UNF4 (4-byte integer)
    // Read as raw bytes
    std::vector<unsigned char> raw_buffer(n_FileSize);
    file.read(reinterpret_cast<char*>(&raw_buffer[0]), n_FileSize);

    IntegerVector result(n_DataSize);
    for (int i = 0; i < n_DataSize; i++) {
      // Manually construct the int32_t with correct byte order (big endian)
      int32_t value = (raw_buffer[i*4] << 24) |
        (raw_buffer[i*4 + 1] << 16) |
        (raw_buffer[i*4 + 2] << 8) |
        raw_buffer[i*4 + 3];
      result[i] = value;
    }
    return result;
  }

  // Default fallback
  return R_NilValue;
}


//'  write_UNF - Write data to UNF file
//' @rdname unf
//' @param data_Export A numeric or integer vector, data to be exported
//' @param fn_UNF A string, name of the UNF file
//' @export
// [[Rcpp::export]]
void write_unf(SEXP data_Export, std::string fn_UNF) {
  // Extract mark from filename
  int mark_UNF = fn_UNF[fn_UNF.length() - 1] - '0';
  int n_Byte = mark_UNF > 0 ? mark_UNF : 4;

  // Open file for writing
  std::ofstream file(fn_UNF, std::ios::binary);
  if (!file.is_open()) {
    stop("Cannot open file for writing: " + fn_UNF);
  }

  // Use n_Byte for type decisions
  if (n_Byte == 4 && mark_UNF == 0) {
    // For UNF0 (4-byte numeric/float)
    NumericVector data = as<NumericVector>(data_Export);
    std::vector<unsigned char> buffer(data.size() * 4);

    for (int i = 0; i < data.size(); i++) {
      // Convert to float and handle endianness
      float value = static_cast<float>(data[i]);
      unsigned char* bytes = reinterpret_cast<unsigned char*>(&value);

      // Write in big-endian order
      buffer[i*4] = bytes[3];
      buffer[i*4 + 1] = bytes[2];
      buffer[i*4 + 2] = bytes[1];
      buffer[i*4 + 3] = bytes[0];
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size());
  }
  else if (n_Byte == 8) {
    // For UNF8 (8-byte numeric/double)
    NumericVector data = as<NumericVector>(data_Export);
    std::vector<unsigned char> buffer(data.size() * 8);

    for (int i = 0; i < data.size(); i++) {
      // Convert to double and handle endianness
      double value = data[i];
      unsigned char* bytes = reinterpret_cast<unsigned char*>(&value);

      // Write in big-endian order
      for (int j = 0; j < 8; j++) {
        buffer[i*8 + j] = bytes[7-j];
      }
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size());
  }
  else if (n_Byte == 1) {
    // For UNF1 (1-byte integer - no swap needed)
    IntegerVector data = as<IntegerVector>(data_Export);
    std::vector<int8_t> buffer(data.size());

    for (int i = 0; i < data.size(); i++) {
      buffer[i] = static_cast<int8_t>(data[i]);
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size());
  }
  else if (n_Byte == 2) {
    // For UNF2 (2-byte integer)
    IntegerVector data = as<IntegerVector>(data_Export);
    std::vector<unsigned char> buffer(data.size() * 2);

    for (int i = 0; i < data.size(); i++) {
      // Get value and convert to big-endian
      int16_t value = static_cast<int16_t>(data[i]);

      // Write bytes in big-endian order
      buffer[i*2] = (value >> 8) & 0xFF;
      buffer[i*2 + 1] = value & 0xFF;
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size());
  }
  else if (n_Byte == 4 && mark_UNF == 4) {
    // For UNF4 (4-byte integer)
    IntegerVector data = as<IntegerVector>(data_Export);
    std::vector<unsigned char> buffer(data.size() * 4);

    for (int i = 0; i < data.size(); i++) {
      // Get value and convert to big-endian
      int32_t value = data[i];

      // Write bytes in big-endian order
      buffer[i*4] = (value >> 24) & 0xFF;
      buffer[i*4 + 1] = (value >> 16) & 0xFF;
      buffer[i*4 + 2] = (value >> 8) & 0xFF;
      buffer[i*4 + 3] = value & 0xFF;
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size());
  }

  file.close();
}
