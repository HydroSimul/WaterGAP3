#include <Rcpp.h>
#include <netcdf.h>
#include <string>
#include <vector>
#include <ctime>
#include <stdexcept>
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// RAII wrapper for NetCDF file closing
struct NCFileCloser {
int ncid;
NCFileCloser(int id) : ncid(id) {}
~NCFileCloser() { nc_close(ncid); }
};

// Function to get current date in YYYY-MM-DD format
std::string getCurrentDate() {
time_t now = time(0);
struct tm tstruct;
char buf[80];
tstruct = *localtime(&now);
strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
return std::string(buf);
}

const std::map<std::string, std::pair<std::string, std::string>> VAR_ATTRIBUTES = {
  {"precipatation_mm", {"mm per m2 per day", "Precipitation"}},
  {"discharge_m3s", {"m3 per second", "Discharge"}},
  {"airTemprature_Cel", {"degree Celsius", "Air Temperature"}},
  {"solarRadition_MJ", {"MJ per m2 per day", "Solar Radiation"}},
  {"netRadiat_MJ", {"MJ per m2 per day", "Net Radiation"}},
  {"relativeHumidity_1", {"ratio (0-1)", "Relative Humidity"}},
  {"vaporPress_kPa", {"kilopascal", "Vapor Pressure"}},
  {"saturatVaporPress_kPa", {"kilopascal", "Saturation Vapor Pressure"}},
  {"windSpeed_m_s", {"meter per second", "Wind Speed"}},
  {"potentialEvatrans_mm", {"mm per m2 per day", "Potential Evapotranspiration"}},
  {"actualEvatrans_mm", {"mm per m2 per day", "Actual Evapotranspiration"}},
  {"runoff_mm", {"mm per m2 per day", "Runoff"}},
  {"rainFall_mm", {"mm per m2 per day", "Precipitation of rainfall"}},
  {"snowFall_mm", {"mm per m2 per day", "Snow Fall"}},
  {"evatrans_mm", {"mm per m2 per day", "Evapotranspiration"}},
  {"soilwater_mm", {"mm  per m2", "Soil Water"}},
  {"groundwater_mm", {"mm  per m2", "Ground Water"}},
  {"soilwater_m3", {"m3", "Soil Water"}},
  {"groundwater_m3", {"m3", "Ground Water"}},
  {"snowIce_mm", {"mm", "Snow and Ice in the snow pack"}},
  {"baseflow_mm", {"mm per m2 per day", "Baseflow"}},
  {"riverwater_m3", {"cubic meter", "River Water"}},
  {"lakewater_m3", {"cubic meter", "Lake Water"}},
  {"lakeEva_mm", {"mm per m2 per day", "Lake Evaporation"}},
  {"riverlakwater_m3", {"cubic meter", "Riverlake Water"}},
  {"riverlakEva_mm", {"mm per m2 per day", "Riverlake Evaporation"}}
};


//' Write data to a NetCDF file in WaterGAP3 format
//'
//' @name nc
//' @param mat_Data_WG3 A numeric matrix containing the data to be written. Rows represent time steps, columns represent spatial units.
//' @param dim_Time An integer vector specifying the time dimension values.
//' @param dim_Spat An integer vector specifying the spatial dimension values.
//' @param path_File A string specifying the directory path where the NetCDF file will be created, for [write_nc_WG3] it is the folder path, for [read_nc_WG3] it should be the file path.
//' @param name_Variable A string specifying the variable name (must be one of: "airTemprature_Cel", "solarRadition_MJ", "precipatation_mm", "discharge_m3s", "runoff_mm").
//' @param str_Continent A string specifying the continent code (must be one of: "eu", "af", "as", "au", "na", "sa").
//' @param suffix_File A string to be appended to the filename (typically indicating time period or version).
//' @return None (creates a NetCDF file on disk) / Matrix
//' @examples
//' \dontrun{
//' data_matrix <- matrix(rnorm(100), nrow=10, ncol=10)
//' write_nc_WG3(data_matrix, 1:10, 1:10, "/path/to/save", "airTemprature_Cel", "eu", "v1")
//' }
//' @export
// [[Rcpp::export]]
void write_nc_WG3(NumericMatrix mat_Data_WG3, IntegerVector dim_Time, IntegerVector dim_Spat,
                 std::string path_File, std::string name_Variable, std::string str_Continent, std::string suffix_File) {
  if (mat_Data_WG3.nrow() != dim_Time.size()) {
    stop("Number of rows in matrix does not match length of time dimension vector.");
  }
  if (mat_Data_WG3.ncol() != dim_Spat.size()) {
    stop("Number of columns in matrix does not match length of spatial dimension vector.");
  }

  std::vector<std::string> valid_continents = {"eu", "af", "as", "au", "na", "sa"};
  if (std::find(valid_continents.begin(), valid_continents.end(), str_Continent) == valid_continents.end()) {
    stop("Invalid continent specified. Must be one of: eu, af, as, au, na, sa");
  }

  std::string filename = path_File + "/" + name_Variable + "_" + str_Continent + "_" + suffix_File + ".nc";

  int ncid, status;
  status = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to create NetCDF file: " + std::string(nc_strerror(status)));
  }

  NCFileCloser fileCloser(ncid);

  int timeDimID, spatDimID;
  status = nc_def_dim(ncid, "time", dim_Time.size(), &timeDimID);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to define time dimension: " + std::string(nc_strerror(status)));
  }
  status = nc_def_dim(ncid, "spat", dim_Spat.size(), &spatDimID);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to define spatial dimension: " + std::string(nc_strerror(status)));
  }

  // Define the time dimension variable
  int timeVarID;
  status = nc_def_var(ncid, "time", NC_INT, 1, &timeDimID, &timeVarID);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to define time dimension variable: " + std::string(nc_strerror(status)));
  }
  // Define the spatial dimension variable
  int spatVarID;
  status = nc_def_var(ncid, "spat", NC_INT, 1, &spatDimID, &spatVarID);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to define spatial dimension variable: " + std::string(nc_strerror(status)));
  }

  int dimids[2] = {timeDimID, spatDimID};
  int dataVarID;
  status = nc_def_var(ncid, name_Variable.c_str(), NC_DOUBLE, 2, dimids, &dataVarID);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to define variable: " + std::string(nc_strerror(status)));
  }

  if (VAR_ATTRIBUTES.find(name_Variable) != VAR_ATTRIBUTES.end()) {
    auto attributes = VAR_ATTRIBUTES.at(name_Variable);
    nc_put_att_text(ncid, dataVarID, "units", attributes.first.size(), attributes.first.c_str());
    nc_put_att_text(ncid, dataVarID, "long_name", attributes.second.size(), attributes.second.c_str());
  }

  nc_put_att_text(ncid, NC_GLOBAL, "creation_date", getCurrentDate().size(), getCurrentDate().c_str());
  nc_put_att_text(ncid, NC_GLOBAL, "institution", 23, "Ruhr-University Bochum");
  nc_put_att_text(ncid, NC_GLOBAL, "creator", 10, "WaterGAP3");
  nc_put_att_text(ncid, NC_GLOBAL, "continent", str_Continent.size(), str_Continent.c_str());

  status = nc_enddef(ncid);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to end NetCDF definition mode: " + std::string(nc_strerror(status)));
  }

  status = nc_put_var_int(ncid, timeVarID, dim_Time.begin());
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to write time dimension variable: " + std::string(nc_strerror(status)));
  }
  status = nc_put_var_int(ncid, spatVarID, dim_Spat.begin());
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to write spatial dimension variable: " + std::string(nc_strerror(status)));
  }

  std::vector<double> data_vector(mat_Data_WG3.nrow() * mat_Data_WG3.ncol());
  for (int i = 0; i < mat_Data_WG3.nrow(); i++) {
    for (int j = 0; j < mat_Data_WG3.ncol(); j++) {
      data_vector[i * mat_Data_WG3.ncol() + j] = mat_Data_WG3(i, j);
    }
  }

  status = nc_put_var_double(ncid, dataVarID, data_vector.data());
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to write data variable: " + std::string(nc_strerror(status)));
  }

  Rcout << "NetCDF file successfully created: " << filename << std::endl;
}

//' Read data from a NetCDF file in WaterGAP3 format
//'
//' @rdname nc
//' @return A numeric matrix containing the data from the NetCDF file, with rows representing time steps and columns representing spatial units.
//' @examples
//' \dontrun{
//' data <- read_nc_WG3("/path/to/file.nc", "airTemprature_Cel")
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix read_nc_WG3(std::string path_File, std::string name_Variable) {
 int ncid, timeDimID, spatDimID, dataVarID;
 int timeLen, spatLen;
 int status = nc_open(path_File.c_str(), NC_NOWRITE, &ncid);
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to open NetCDF file: " + std::string(nc_strerror(status)));
 }

 NCFileCloser fileCloser(ncid);

 status = nc_inq_dimid(ncid, "time", &timeDimID);
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to get time dimension ID: " + std::string(nc_strerror(status)));
 }
 status = nc_inq_dimlen(ncid, timeDimID, (size_t*)&timeLen);
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to get time dimension length: " + std::string(nc_strerror(status)));
 }

 status = nc_inq_dimid(ncid, "spat", &spatDimID);
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to get spatial dimension ID: " + std::string(nc_strerror(status)));
 }
 status = nc_inq_dimlen(ncid, spatDimID, (size_t*)&spatLen);
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to get spatial dimension length: " + std::string(nc_strerror(status)));
 }

 status = nc_inq_varid(ncid, name_Variable.c_str(), &dataVarID);
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to get variable ID for " + name_Variable + ": " + std::string(nc_strerror(status)));
 }

 std::vector<double> data_vector(timeLen * spatLen);
 status = nc_get_var_double(ncid, dataVarID, data_vector.data());
 if (status != NC_NOERR) {
   throw std::runtime_error("Failed to read variable data: " + std::string(nc_strerror(status)));
 }

 NumericMatrix mat_Data_WG3(timeLen, spatLen);
 for (int i = 0; i < timeLen; i++) {
   for (int j = 0; j < spatLen; j++) {
     mat_Data_WG3(i, j) = data_vector[i * spatLen + j];
   }
 }

 return mat_Data_WG3;
}


//' @rdname nc
//' @return A numeric vector containing the data from the NetCDF file, with rows representing time steps and columns representing spatial units.
//' @examples
//' \dontrun{
//' data <- read_nc_dim_WG3("/path/to/file.nc", "time")
//' }
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector read_nc_dim_WG3(std::string path_File, std::string dim_name) {
  int ncid, dimID, varID;
  size_t dimLen;
  int status = nc_open(path_File.c_str(), NC_NOWRITE, &ncid);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to open NetCDF file: " + std::string(nc_strerror(status)));
  }
  NCFileCloser fileCloser(ncid);

  // Get dimension ID and length
  status = nc_inq_dimid(ncid, dim_name.c_str(), &dimID);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to get dimension ID for " + dim_name + ": " +
                             std::string(nc_strerror(status)));
  }

  status = nc_inq_dimlen(ncid, dimID, &dimLen);
  if (status != NC_NOERR) {
    throw std::runtime_error("Failed to get dimension length for " + dim_name + ": " +
                             std::string(nc_strerror(status)));
  }

  // Create an IntegerVector with indices
  Rcpp::IntegerVector dim_indices(dimLen);

  // Try to get dimension variable values if it exists
  status = nc_inq_varid(ncid, dim_name.c_str(), &varID);
  if (status == NC_NOERR) {
    // Dimension variable exists, read it
    std::vector<int> dim_values(dimLen);
    status = nc_get_var_int(ncid, varID, dim_values.data());
    if (status == NC_NOERR) {
      // Copy values to the IntegerVector
      for (size_t i = 0; i < dimLen; i++) {
        dim_indices[i] = dim_values[i];
      }
    } else {
      // Try reading as doubles and converting to int if int read fails
      std::vector<double> double_values(dimLen);
      status = nc_get_var_double(ncid, varID, double_values.data());
      if (status == NC_NOERR) {
        for (size_t i = 0; i < dimLen; i++) {
          dim_indices[i] = static_cast<int>(double_values[i]);
        }
      } else {
        // If both fail, use sequence indices
        for (size_t i = 0; i < dimLen; i++) {
          dim_indices[i] = i;
        }
      }
    }
  } else {
    // No dimension variable exists, use sequence indices
    for (size_t i = 0; i < dimLen; i++) {
      dim_indices[i] = i;
    }
  }

  return dim_indices;
}


//' @rdname nc
//'
//' @description Combine multiple NetCDF files along the time dimension
//'
//' This function reads multiple NetCDF files with the same spatial dimension but different
//' time steps, combines them along the time dimension, and writes the result to a new NetCDF file.
//'
//' @param path_FilesBind A character vector of paths to the input NetCDF files
//' @param path_Out Path to the output NetCDF file
//' @param name_Variable Name of the variable to combine
//' @return None (writes a NetCDF file to disk)
//' @export
// [[Rcpp::export]]
void bind_nc_WG3(std::vector<std::string> path_FilesBind,
                std::string path_Out,
                std::string name_Variable) {
 if (path_FilesBind.empty()) {
   stop("No input files provided");
 }

 // Open first file to get dimensions and attributes
 int ncid_first, status;
 status = nc_open(path_FilesBind[0].c_str(), NC_NOWRITE, &ncid_first);
 if (status != NC_NOERR) {
   stop("Failed to open first NetCDF file: " + std::string(nc_strerror(status)));
 }
 NCFileCloser firstCloser(ncid_first);

 // Get dimensions from first file
 int timeDimID_first, spatDimID_first;
 size_t timeLen_first, spatLen_first;

 status = nc_inq_dimid(ncid_first, "time", &timeDimID_first);
 if (status != NC_NOERR) {
   stop("Failed to get time dimension ID from first file: " + std::string(nc_strerror(status)));
 }
 status = nc_inq_dimlen(ncid_first, timeDimID_first, &timeLen_first);
 if (status != NC_NOERR) {
   stop("Failed to get time dimension length from first file: " + std::string(nc_strerror(status)));
 }

 status = nc_inq_dimid(ncid_first, "spat", &spatDimID_first);
 if (status != NC_NOERR) {
   stop("Failed to get spatial dimension ID from first file: " + std::string(nc_strerror(status)));
 }
 status = nc_inq_dimlen(ncid_first, spatDimID_first, &spatLen_first);
 if (status != NC_NOERR) {
   stop("Failed to get spatial dimension length from first file: " + std::string(nc_strerror(status)));
 }

 // Get time and spatial values from first file
 std::vector<int> all_time_values;
 std::vector<int> time_values_first(timeLen_first);
 int timeVarID_first;
 status = nc_inq_varid(ncid_first, "time", &timeVarID_first);
 if (status == NC_NOERR) {
   status = nc_get_var_int(ncid_first, timeVarID_first, time_values_first.data());
   if (status != NC_NOERR) {
     stop("Failed to read time values from first file: " + std::string(nc_strerror(status)));
   }
   all_time_values.insert(all_time_values.end(), time_values_first.begin(), time_values_first.end());
 } else {
   // If time variable doesn't exist, create sequential values
   for (size_t i = 0; i < timeLen_first; i++) {
     all_time_values.push_back(i);
   }
 }

 std::vector<int> spat_values(spatLen_first);
 int spatVarID_first;
 status = nc_inq_varid(ncid_first, "spat", &spatVarID_first);
 if (status == NC_NOERR) {
   status = nc_get_var_int(ncid_first, spatVarID_first, spat_values.data());
   if (status != NC_NOERR) {
     stop("Failed to read spatial values from first file: " + std::string(nc_strerror(status)));
   }
 } else {
   // If spatial variable doesn't exist, create sequential values
   for (size_t i = 0; i < spatLen_first; i++) {
     spat_values[i] = i;
   }
 }

 // Read data from first file
 int dataVarID_first;
 status = nc_inq_varid(ncid_first, name_Variable.c_str(), &dataVarID_first);
 if (status != NC_NOERR) {
   stop("Failed to get variable ID from first file: " + std::string(nc_strerror(status)));
 }

 std::vector<double> combined_data(timeLen_first * spatLen_first);
 status = nc_get_var_double(ncid_first, dataVarID_first, combined_data.data());
 if (status != NC_NOERR) {
   stop("Failed to read data from first file: " + std::string(nc_strerror(status)));
 }

 // Process remaining files
 for (size_t f = 1; f < path_FilesBind.size(); f++) {
   int ncid, status;
   status = nc_open(path_FilesBind[f].c_str(), NC_NOWRITE, &ncid);
   if (status != NC_NOERR) {
     stop("Failed to open NetCDF file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }
   NCFileCloser closer(ncid);

   // Check dimensions match
   int timeDimID, spatDimID;
   size_t timeLen, spatLen;

   status = nc_inq_dimid(ncid, "time", &timeDimID);
   if (status != NC_NOERR) {
     stop("Failed to get time dimension ID from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }
   status = nc_inq_dimlen(ncid, timeDimID, &timeLen);
   if (status != NC_NOERR) {
     stop("Failed to get time dimension length from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }

   status = nc_inq_dimid(ncid, "spat", &spatDimID);
   if (status != NC_NOERR) {
     stop("Failed to get spatial dimension ID from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }
   status = nc_inq_dimlen(ncid, spatDimID, &spatLen);
   if (status != NC_NOERR) {
     stop("Failed to get spatial dimension length from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }

   if (spatLen != spatLen_first) {
     stop("Spatial dimension length mismatch in file " + path_FilesBind[f]);
   }

   // Get time values
   std::vector<int> time_values(timeLen);
   int timeVarID;
   status = nc_inq_varid(ncid, "time", &timeVarID);
   if (status == NC_NOERR) {
     status = nc_get_var_int(ncid, timeVarID, time_values.data());
     if (status != NC_NOERR) {
       stop("Failed to read time values from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
     }
     all_time_values.insert(all_time_values.end(), time_values.begin(), time_values.end());
   } else {
     // If time variable doesn't exist, create sequential values
     size_t current_size = all_time_values.size();
     for (size_t i = 0; i < timeLen; i++) {
       all_time_values.push_back(current_size + i);
     }
   }

   // Read data
   int dataVarID;
   status = nc_inq_varid(ncid, name_Variable.c_str(), &dataVarID);
   if (status != NC_NOERR) {
     stop("Failed to get variable ID from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }

   std::vector<double> file_data(timeLen * spatLen);
   status = nc_get_var_double(ncid, dataVarID, file_data.data());
   if (status != NC_NOERR) {
     stop("Failed to read data from file " + path_FilesBind[f] + ": " + std::string(nc_strerror(status)));
   }

   // Append data to combined_data
   combined_data.insert(combined_data.end(), file_data.begin(), file_data.end());
 }

 // Create output file
 int ncid_out, status_out;
 status_out = nc_create(path_Out.c_str(), NC_CLOBBER, &ncid_out);
 if (status_out != NC_NOERR) {
   stop("Failed to create output NetCDF file: " + std::string(nc_strerror(status_out)));
 }
 NCFileCloser outCloser(ncid_out);

 // Define dimensions
 int timeDimID_out, spatDimID_out;
 size_t total_time_len = all_time_values.size();

 status_out = nc_def_dim(ncid_out, "time", total_time_len, &timeDimID_out);
 if (status_out != NC_NOERR) {
   stop("Failed to define time dimension in output file: " + std::string(nc_strerror(status_out)));
 }

 status_out = nc_def_dim(ncid_out, "spat", spatLen_first, &spatDimID_out);
 if (status_out != NC_NOERR) {
   stop("Failed to define spatial dimension in output file: " + std::string(nc_strerror(status_out)));
 }

 // Define variables
 int timeVarID_out;
 status_out = nc_def_var(ncid_out, "time", NC_INT, 1, &timeDimID_out, &timeVarID_out);
 if (status_out != NC_NOERR) {
   stop("Failed to define time variable in output file: " + std::string(nc_strerror(status_out)));
 }

 int spatVarID_out;
 status_out = nc_def_var(ncid_out, "spat", NC_INT, 1, &spatDimID_out, &spatVarID_out);
 if (status_out != NC_NOERR) {
   stop("Failed to define spatial variable in output file: " + std::string(nc_strerror(status_out)));
 }

 int dimids[2] = {timeDimID_out, spatDimID_out};
 int dataVarID_out;
 status_out = nc_def_var(ncid_out, name_Variable.c_str(), NC_DOUBLE, 2, dimids, &dataVarID_out);
 if (status_out != NC_NOERR) {
   stop("Failed to define data variable in output file: " + std::string(nc_strerror(status_out)));
 }

 // Copy attributes from first file's data variable
 int dataVarID_first_attrs;
 status = nc_inq_varid(ncid_first, name_Variable.c_str(), &dataVarID_first_attrs);
 if (status == NC_NOERR) {
   // Copy units attribute if it exists
   size_t units_len;
   status = nc_inq_attlen(ncid_first, dataVarID_first_attrs, "units", &units_len);
   if (status == NC_NOERR) {
     std::vector<char> units(units_len + 1);
     status = nc_get_att_text(ncid_first, dataVarID_first_attrs, "units", units.data());
     if (status == NC_NOERR) {
       units[units_len] = '\0';
       status_out = nc_put_att_text(ncid_out, dataVarID_out, "units", units_len, units.data());
       if (status_out != NC_NOERR) {
         Rcerr << "Warning: Failed to copy units attribute to output file: " << nc_strerror(status_out) << std::endl;
       }
     }
   }

   // Copy long_name attribute if it exists
   size_t long_name_len;
   status = nc_inq_attlen(ncid_first, dataVarID_first_attrs, "long_name", &long_name_len);
   if (status == NC_NOERR) {
     std::vector<char> long_name(long_name_len + 1);
     status = nc_get_att_text(ncid_first, dataVarID_first_attrs, "long_name", long_name.data());
     if (status == NC_NOERR) {
       long_name[long_name_len] = '\0';
       status_out = nc_put_att_text(ncid_out, dataVarID_out, "long_name", long_name_len, long_name.data());
       if (status_out != NC_NOERR) {
         Rcerr << "Warning: Failed to copy long_name attribute to output file: " << nc_strerror(status_out) << std::endl;
       }
     }
   }
 }

 // Copy global attributes from first file
 int ngatts;
 status = nc_inq_natts(ncid_first, &ngatts);
 if (status == NC_NOERR) {
   for (int i = 0; i < ngatts; i++) {
     char name[NC_MAX_NAME + 1];
     status = nc_inq_attname(ncid_first, NC_GLOBAL, i, name);
     if (status == NC_NOERR) {
       nc_type type;
       size_t len;
       status = nc_inq_att(ncid_first, NC_GLOBAL, name, &type, &len);
       if (status == NC_NOERR) {
         if (type == NC_CHAR) {
           std::vector<char> value(len + 1);
           status = nc_get_att_text(ncid_first, NC_GLOBAL, name, value.data());
           if (status == NC_NOERR) {
             value[len] = '\0';
             status_out = nc_put_att_text(ncid_out, NC_GLOBAL, name, len, value.data());
           }
         } else if (type == NC_INT) {
           std::vector<int> value(len);
           status = nc_get_att_int(ncid_first, NC_GLOBAL, name, value.data());
           if (status == NC_NOERR) {
             status_out = nc_put_att_int(ncid_out, NC_GLOBAL, name, type, len, value.data());
           }
         } else if (type == NC_DOUBLE) {
           std::vector<double> value(len);
           status = nc_get_att_double(ncid_first, NC_GLOBAL, name, value.data());
           if (status == NC_NOERR) {
             status_out = nc_put_att_double(ncid_out, NC_GLOBAL, name, type, len, value.data());
           }
         }

         if (status_out != NC_NOERR) {
           Rcerr << "Warning: Failed to copy global attribute " << name << " to output file: " << nc_strerror(status_out) << std::endl;
         }
       }
     }
   }
 }

 status_out = nc_enddef(ncid_out);
 if (status_out != NC_NOERR) {
   stop("Failed to end definition mode in output file: " + std::string(nc_strerror(status_out)));
 }

 // Write data to output file
 status_out = nc_put_var_int(ncid_out, timeVarID_out, all_time_values.data());
 if (status_out != NC_NOERR) {
   stop("Failed to write time values to output file: " + std::string(nc_strerror(status_out)));
 }

 status_out = nc_put_var_int(ncid_out, spatVarID_out, spat_values.data());
 if (status_out != NC_NOERR) {
   stop("Failed to write spatial values to output file: " + std::string(nc_strerror(status_out)));
 }

 status_out = nc_put_var_double(ncid_out, dataVarID_out, combined_data.data());
 if (status_out != NC_NOERR) {
   stop("Failed to write data to output file: " + std::string(nc_strerror(status_out)));
 }

 Rcout << "Successfully combined " << path_FilesBind.size() << " files into " << path_Out << std::endl;
}





