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

  // Function to swap endianness
  auto swap_endian = [](unsigned int x) {
    return ((x & 0x000000FF) << 24) |
      ((x & 0x0000FF00) << 8) |
      ((x & 0x00FF0000) >> 8) |
      ((x & 0xFF000000) >> 24);
  };

  auto swap_endian_short = [](unsigned short x) {
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
        result[i] = buffer[i];
      }
      return result;
    }
    else if (n_Byte == 2) {
      // For UNF2 (2-byte integer)
      std::vector<int16_t> buffer(n_DataSize);
      file.read(reinterpret_cast<char*>(&buffer[0]), n_FileSize);

      IntegerVector result(n_DataSize);
      for (int i = 0; i < n_DataSize; i++) {
        result[i] = swap_endian_short(buffer[i]);
      }
      return result;
    }
  }
  else if (type_UNF == "integer") {
    // For UNF4 (4-byte integer)
    std::vector<int32_t> buffer(n_DataSize);
    file.read(reinterpret_cast<char*>(&buffer[0]), n_FileSize);

    IntegerVector result(n_DataSize);
    for (int i = 0; i < n_DataSize; i++) {
      result[i] = swap_endian(buffer[i]);
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

  // Endian swap functions
  auto swap_endian = [](unsigned int x) {
    return ((x & 0x000000FF) << 24) |
      ((x & 0x0000FF00) << 8) |
      ((x & 0x00FF0000) >> 8) |
      ((x & 0xFF000000) >> 24);
  };

  auto swap_endian_short = [](unsigned short x) {
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

  // Use n_Byte for type decisions
  if (n_Byte == 4 && mark_UNF == 0) {
    // For UNF0 (4-byte numeric/float)
    NumericVector data = as<NumericVector>(data_Export);
    std::vector<float> buffer(data.size());

    for (int i = 0; i < data.size(); i++) {
      buffer[i] = swap_endian_float(static_cast<float>(data[i]));
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size() * sizeof(float));
  }
  else if (n_Byte == 8) {
    // For UNF8 (8-byte numeric/double)
    NumericVector data = as<NumericVector>(data_Export);
    std::vector<double> buffer(data.size());

    for (int i = 0; i < data.size(); i++) {
      buffer[i] = swap_endian_double(data[i]);
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size() * sizeof(double));
  }
  else if (n_Byte == 1) {
    // For UNF1 (1-byte integer - no swap needed)
    IntegerVector data = as<IntegerVector>(data_Export);
    std::vector<int8_t> buffer(data.size());

    for (int i = 0; i < data.size(); i++) {
      buffer[i] = static_cast<int8_t>(data[i]);
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size() * sizeof(int8_t));
  }
  else if (n_Byte == 2) {
    // For UNF2 (2-byte integer)
    IntegerVector data = as<IntegerVector>(data_Export);
    std::vector<int16_t> buffer(data.size());

    for (int i = 0; i < data.size(); i++) {
      buffer[i] = swap_endian_short(static_cast<int16_t>(data[i]));
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size() * sizeof(int16_t));
  }
  else if (n_Byte == 4 && mark_UNF == 4) {
    // For UNF4 (4-byte integer)
    IntegerVector data = as<IntegerVector>(data_Export);
    std::vector<int32_t> buffer(data.size());

    for (int i = 0; i < data.size(); i++) {
      buffer[i] = swap_endian(data[i]);
    }

    file.write(reinterpret_cast<char*>(&buffer[0]), buffer.size() * sizeof(int32_t));
  }

  file.close();
}
