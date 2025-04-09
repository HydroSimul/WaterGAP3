#include "model.h"
// [[Rcpp::interfaces(r, cpp)]]

//' run_WaterGAP3
//' @name run_WaterGAP3
//' @description
//' Run WaterGAP3_N hydrological model
//'
//' This function runs the WaterGAP3_N hydrological model for a specified region
//' and time period, returning discharge calculations.
//'
//' @inheritParams WaterGAP3
//' @param path_MeteoInput Directory path containing meteorological input files.
//' @param path_HydroParam Directory path containing hydrological parameter files.
//' @param path_InitialState Directory path containing initial state files. If "UNKNOW", default initial states will be used.
//' @param path_Boundary Directory path containing boundary condition files. If "UNKNOW", default boundary conditions will be used.
//'
//' @return A numeric matrix (CELL_discharge_m3) containing discharge values for each cell and time step.
//' @export
// [[Rcpp::export]]
NumericMatrix run_WaterGAP3_N(std::string name_Region, std::string mark_Time, int n_Time, int n_Spat,
                                  std::string path_MeteoInput, std::string path_HydroParam,
                                  std::string path_InitialState = "UNKNOW",
                                  std::string path_Boundary = "UNKNOW",
                                  std::string path_VariExport = "NonExport",
                                  std::string path_FinalState = "NonExport") {

 // Meteo Forcing Matrix
 NumericMatrix ATMOS_precipitation_mm = load_wgmat(path_MeteoInput + "ATMOS_precipitation_mm_" + name_Region + "_" + mark_Time + ".wgmat");
 NumericMatrix ATMOS_temperature_Cel = load_wgmat(path_MeteoInput + "ATMOS_temperature_Cel_" + name_Region + "_" + mark_Time + ".wgmat");
 NumericMatrix ATMOS_solarRadiat_MJ = load_wgmat(path_MeteoInput + "ATMOS_solarRadiat_MJ_" + name_Region + "_" + mark_Time + ".wgmat");
 NumericMatrix ATMOS_solarRadiatClearSky_MJ = load_wgmat(path_MeteoInput + "ATMOS_solarRadiatClearSky_MJ_" + name_Region + "_" + mark_Time + ".wgmat");

 // Calculate LAND_interceptCapacity_mm
 NumericMatrix LAND_leafAreaIndex_1 = load_wgmat(path_MeteoInput + "LAND_leafAreaIndex_1_" + name_Region + "_" + mark_Time + ".wgmat");
 NumericMatrix LAND_interceptCapacity_mm = LAND_leafAreaIndex_1 * 0.3 + 0.01;

 // Upstream boundary

 // Upstream boundary
 IntegerVector Upstream_cellNumber_int;
 NumericMatrix Upstream_streamflow_m3;

 if (path_Boundary == "UNKNOW") {
   Upstream_cellNumber_int = IntegerVector::create(0);
   // Use NumericMatrix constructor properly
   Upstream_streamflow_m3 = NumericMatrix(1, 1);
   Upstream_streamflow_m3(0, 0) = 0.0;
 } else {
   Upstream_cellNumber_int = read_unf(path_Boundary + "Riverlak_cellNumber_int_" + name_Region + ".UNF4");
   Upstream_streamflow_m3 = load_wgmat(path_Boundary + "Upstream_streamflow_m3_" + name_Region + ".wgmat");
 }
 // Hydro Parameter Vector
 // LAND
 NumericVector LAND_area_km2 = read_unf(path_HydroParam + "LAND_area_km2_" + name_Region + ".UNF0");
 NumericVector LAND_albedo_1 = read_unf(path_HydroParam + "LAND_albedo_1_" + name_Region + ".UNF0");
 NumericVector LAND_snowAlbedo_1 = read_unf(path_HydroParam + "LAND_snowAlbedo_1_" + name_Region + ".UNF0");
 NumericVector LAND_builtRatio_1 = read_unf(path_HydroParam + "LAND_builtRatio_1_" + name_Region + ".UNF0");

 // SOIL
 NumericVector SOIL_capacity_mm = read_unf(path_HydroParam + "SOIL_capacity_mm_" + name_Region + ".UNF0");
 NumericVector SOIL_potentialPercola_mm = read_unf(path_HydroParam + "SOIL_potentialPercola_mm_" + name_Region + ".UNF0");

 // RIVER
 NumericVector RIVER_length_km = read_unf(path_HydroParam + "RIVER_length_km_" + name_Region + ".UNF0");
 NumericVector RIVER_velocity_km = read_unf(path_HydroParam + "RIVER_velocity_km_" + name_Region + ".UNF0");

 // CELL
 NumericVector CELL_elevation_m = read_unf(path_HydroParam + "CELL_elevation_m_" + name_Region + ".UNF0");
 List CELL_cellNumberStep_int = read_int_vector_list(path_HydroParam + "CELL_cellNumberStep_int_" + name_Region + ".bin");
 List CELL_inflowCellNumberStep_int = read_int_matrix_list(path_HydroParam + "CELL_inflowCellNumberStep_int_" + name_Region + ".bin");

 // Lake
 IntegerVector Lake_cellNumber_int = read_unf(path_HydroParam + "Lake_cellNumber_int_" + name_Region + ".UNF4");
 NumericVector Lake_area_km2 = read_unf(path_HydroParam + "Lake_area_km2_" + name_Region + ".UNF0");
 NumericVector Lake_capacity_m3 = read_unf(path_HydroParam + "Lake_capacity_m3_" + name_Region + ".UNF0");
 NumericVector Lake_albedo_1 = read_unf(path_HydroParam + "Lake_albedo_1_" + name_Region + ".UNF0");

 // Riverlak
 IntegerVector Riverlak_cellNumber_int = read_unf(path_HydroParam + "Riverlak_cellNumber_int_" + name_Region + ".UNF4");
 NumericVector Riverlak_area_km2 = read_unf(path_HydroParam + "Riverlak_area_km2_" + name_Region + ".UNF0");
 NumericVector Riverlak_capacity_m3 = read_unf(path_HydroParam + "Riverlak_capacity_m3_" + name_Region + ".UNF0");
 NumericVector Riverlak_albedo_1 = read_unf(path_HydroParam + "Riverlak_albedo_1_" + name_Region + ".UNF0");

 // Initial states
 NumericVector SNOW_ice_mm;
 NumericVector LAND_interceptWater_mm;
 NumericVector SOIL_water_mm;
 NumericVector GROUND_water_mm;
 NumericVector RIVER_water_m3;
 NumericVector Lake_water_m3;
 NumericVector Riverlak_water_m3;

 if (path_InitialState == "UNKNOW") {
   SNOW_ice_mm = NumericVector(n_Spat, 0.0);
   LAND_interceptWater_mm = NumericVector(n_Spat, 0.0);
   SOIL_water_mm = SOIL_capacity_mm * 0.6;
   GROUND_water_mm = NumericVector(n_Spat, 0.1);
   RIVER_water_m3 = NumericVector(n_Spat, 0.1);
   Lake_water_m3 = Lake_capacity_m3 * 0.6;
   Riverlak_water_m3 = Riverlak_capacity_m3 * 0.6;
 } else {
   SNOW_ice_mm = read_unf(path_InitialState + "SNOW_ice_mm_" + name_Region + ".UNF0");
   LAND_interceptWater_mm = read_unf(path_InitialState + "LAND_interceptWater_mm_" + name_Region + ".UNF0");
   SOIL_water_mm = read_unf(path_InitialState + "SOIL_water_mm_" + name_Region + ".UNF0");
   GROUND_water_mm = read_unf(path_InitialState + "GROUND_water_mm_" + name_Region + ".UNF0");
   RIVER_water_m3 = read_unf(path_InitialState + "RIVER_water_m3_" + name_Region + ".UNF0");
   Lake_water_m3 = read_unf(path_InitialState + "Lake_water_m3_" + name_Region + ".UNF0");
   Riverlak_water_m3 = read_unf(path_InitialState + "Riverlak_water_m3_" + name_Region + ".UNF0");
 }

 // Parameters
 NumericVector param_ATMOS_thr_Ts = read_unf(path_HydroParam + "param_ATMOS_thr_Ts_" + name_Region + ".UNF0");
 NumericVector param_SNOW_fac_f = read_unf(path_HydroParam + "param_SNOW_fac_f_" + name_Region + ".UNF0");
 NumericVector param_SNOW_fac_Tmelt = read_unf(path_HydroParam + "param_SNOW_fac_Tmelt_" + name_Region + ".UNF0");
 NumericVector param_EVATRANS_prt_alpha = read_unf(path_HydroParam + "param_EVATRANS_prt_alpha_" + name_Region + ".UNF0");
 NumericVector param_EVATRANS_vic_gamma = read_unf(path_HydroParam + "param_EVATRANS_vic_gamma_" + name_Region + ".UNF0");
 NumericVector param_EVATRANS_sup_k = read_unf(path_HydroParam + "param_EVATRANS_sup_k_" + name_Region + ".UNF0");
 NumericVector param_EVATRANS_sup_gamma = read_unf(path_HydroParam + "param_EVATRANS_sup_gamma_" + name_Region + ".UNF0");
 NumericVector param_EVATRANS_wat_petmax = read_unf(path_HydroParam + "param_EVATRANS_wat_petmax_" + name_Region + ".UNF0");
 NumericVector param_INFILT_hbv_beta = read_unf(path_HydroParam + "param_INFILT_hbv_beta_" + name_Region + ".UNF0");
 LogicalVector param_PERCOLA_wat_01 = read_unf(path_HydroParam + "param_PERCOLA_wat_01_" + name_Region + ".UNF1");
 NumericVector param_PERCOLA_wat_k = read_unf(path_HydroParam + "param_PERCOLA_wat_k_" + name_Region + ".UNF0");
 NumericVector param_PERCOLA_wat_thresh = read_unf(path_HydroParam + "param_PERCOLA_wat_thresh_" + name_Region + ".UNF0");
 NumericVector param_BASEFLOW_sur_k = read_unf(path_HydroParam + "param_BASEFLOW_sur_k_" + name_Region + ".UNF0");
 NumericVector param_Lake_acp_storeFactor = read_unf(path_HydroParam + "param_Lake_acp_storeFactor_" + name_Region + ".UNF0");
 NumericVector param_Lake_acp_gamma = read_unf(path_HydroParam + "param_Lake_acp_gamma_" + name_Region + ".UNF0");
 NumericVector param_Riverlak_lin_storeFactor = read_unf(path_HydroParam + "param_Riverlak_lin_storeFactor_" + name_Region + ".UNF0");

 // Call the main WaterGAP3_N function (this would need to be implemented in C++ as well)
 NumericMatrix CELL_discharge_m3 = WaterGAP3_N(name_Region,
                                               mark_Time,
                                               n_Time,
                                               n_Spat,
                                               ATMOS_precipitation_mm,
                                               ATMOS_temperature_Cel,
                                               ATMOS_solarRadiat_MJ,
                                               ATMOS_solarRadiatClearSky_MJ,
                                               Upstream_cellNumber_int,
                                               Upstream_streamflow_m3,
                                               SNOW_ice_mm,
                                               LAND_area_km2,
                                               LAND_albedo_1,
                                               LAND_snowAlbedo_1,
                                               LAND_builtRatio_1,
                                               LAND_interceptWater_mm,
                                               LAND_interceptCapacity_mm,
                                               SOIL_water_mm,
                                               SOIL_capacity_mm,
                                               SOIL_potentialPercola_mm,
                                               GROUND_water_mm,
                                               RIVER_water_m3,
                                               RIVER_length_km,
                                               RIVER_velocity_km,
                                               CELL_elevation_m,
                                               CELL_cellNumberStep_int,
                                               CELL_inflowCellNumberStep_int,
                                               Lake_cellNumber_int,
                                               Lake_water_m3,
                                               Lake_area_km2,
                                               Lake_capacity_m3,
                                               Lake_albedo_1,
                                               Riverlak_cellNumber_int,
                                               Riverlak_water_m3,
                                               Riverlak_area_km2,
                                               Riverlak_capacity_m3,
                                               Riverlak_albedo_1,
                                               param_ATMOS_thr_Ts,
                                               param_SNOW_fac_f,
                                               param_SNOW_fac_Tmelt,
                                               param_EVATRANS_prt_alpha,
                                               param_EVATRANS_vic_gamma,
                                               param_EVATRANS_sup_k,
                                               param_EVATRANS_sup_gamma,
                                               param_EVATRANS_wat_petmax,
                                               param_INFILT_hbv_beta,
                                               param_PERCOLA_wat_01,
                                               param_PERCOLA_wat_k,
                                               param_PERCOLA_wat_thresh,
                                               param_BASEFLOW_sur_k,
                                               param_Lake_acp_storeFactor,
                                               param_Lake_acp_gamma,
                                               param_Riverlak_lin_storeFactor,
                                               path_FinalState,
                                               path_VariExport);

 // Return results
 return CELL_discharge_m3;
}
