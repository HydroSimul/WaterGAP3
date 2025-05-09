#include "model.h"
// [[Rcpp::interfaces(r, cpp)]]
//' run_model
//' @name run_model
//' @description
//' Run WaterGAP3_N hydrological model
//'
//' This function runs the WaterGAP3_N hydrological model for a specified region
//' and time period, returning discharge calculations.
//'
//' @inheritParams model
//' @param path_MeteoInput Directory path containing meteorological input files.
//' @param path_WateruseInput Directory path containing wateruse input files.
//' @param path_HydroParam Directory path containing hydrological parameter files.
//' @param path_InitialState Directory path containing initial state files. If "UNKNOW", default initial states will be used.
//' @param path_Boundary Directory path containing boundary condition files. If "UNKNOW", default boundary conditions will be used.
//'
//' @return A numeric matrix (CELL_discharge_m3) containing discharge values for each cell and time step.
//' @export
// [[Rcpp::export]]
NumericMatrix run_WaterGAP3_U(std::string name_Region, std::string mark_Time, int n_Time, int n_Spat,
                           std::string path_MeteoInput, std::string path_HydroParam, std::string path_WateruseInput,
                           std::string path_InitialState = "UNKNOW",
                           std::string path_Boundary = "UNKNOW",
                           std::string path_VariExport = "NonExport",
                           std::string path_FinalState = "NonExport") {


 const std::vector<int> month_days_normal = {31, 28, 31, 30, 31, 30,
                                             31, 31, 30, 31, 30, 31};
 const std::vector<int> month_days_leap = {31, 29, 31, 30, 31, 30,
                                           31, 31, 30, 31, 30, 31};
 const std::vector<int>& days_in_month = (n_Time == 365) ? month_days_normal : month_days_leap;
 std::vector<int> month_vec;
 month_vec.reserve(std::accumulate(days_in_month.begin(), days_in_month.end(), 0));
 for (size_t m = 0; m < days_in_month.size(); ++m) {
   month_vec.insert(month_vec.end(), days_in_month[m], m + 1);
 }
 IntegerVector idx_Month = wrap(month_vec);


 // Meteo Forcing Matrix
 NumericMatrix ATMOS_precipitation_mm = load_matbin(path_MeteoInput + "ATMOS_precipitation_mm" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix ATMOS_temperature_Cel = load_matbin(path_MeteoInput + "ATMOS_temperature_Cel" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix ATMOS_solarRadiat_MJ = load_matbin(path_MeteoInput + "ATMOS_solarRadiat_MJ" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix ATMOS_solarRadiatClearSky_MJ = load_matbin(path_MeteoInput + "ATMOS_solarRadiatClearSky_MJ" + "_" + name_Region + "_" + mark_Time + ".matbin");

 // Withdraw
 NumericMatrix WITHDRAW_surface_m3 = load_matbin(path_WateruseInput + "WITHDRAW_surface_m3" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix WITHDRAW_ground_m3 = load_matbin(path_WateruseInput + "WITHDRAW_ground_m3" + "_" + name_Region + "_" + mark_Time + ".matbin");



 // Calculate LAND_interceptCapacity_mm
 NumericMatrix LAND_leafAreaIndex_1 = load_matbin(path_MeteoInput + "LAND_leafAreaIndex_1" + "_" + name_Region + "_" + mark_Time + ".matbin");
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
   Upstream_cellNumber_int = load_vecbin(path_Boundary + "Upstream_cellNumber_int" + "_" + name_Region + ".vecbin");
   Upstream_streamflow_m3 = load_matbin(path_Boundary + "Upstream_streamflow_m3" + "_" + name_Region + ".matbin");
 }
 // Hydro Parameter Vector
 // LAND
 NumericVector LAND_area_km2 = load_vecbin(path_HydroParam + "LAND_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector LAND_albedo_1 = load_vecbin(path_HydroParam + "LAND_albedo_1" + "_" + name_Region + ".vecbin");
 NumericVector LAND_snowAlbedo_1 = load_vecbin(path_HydroParam + "LAND_snowAlbedo_1" + "_" + name_Region + ".vecbin");
 NumericVector LAND_builtRatio_1 = load_vecbin(path_HydroParam + "LAND_builtRatio_1" + "_" + name_Region + ".vecbin");

 // SOIL
 NumericVector SOIL_capacity_mm = load_vecbin(path_HydroParam + "SOIL_capacity_mm" + "_" + name_Region + ".vecbin");
 NumericVector SOIL_potentialPercola_mm = load_vecbin(path_HydroParam + "SOIL_potentialPercola_mm" + "_" + name_Region + ".vecbin");

 // RIVER
 NumericVector RIVER_length_km = load_vecbin(path_HydroParam + "RIVER_length_km" + "_" + name_Region + ".vecbin");
 NumericVector RIVER_velocity_km = load_vecbin(path_HydroParam + "RIVER_velocity_km" + "_" + name_Region + ".vecbin");

 // CELL
 NumericVector CELL_elevation_m = load_vecbin(path_HydroParam + "CELL_elevation_m" + "_" + name_Region + ".vecbin");
 List CELL_cellNumberStep_int = read_int_vector_list(path_HydroParam + "CELL_cellNumberStep_int" + "_" + name_Region + ".bin");
 List CELL_inflowCellNumberStep_int = read_int_matrix_list(path_HydroParam + "CELL_inflowCellNumberStep_int" + "_" + name_Region + ".bin");
 IntegerMatrix CELL_cellNumberAround_int = load_matbin(path_HydroParam + "CELL_elevation_m" + "_" + name_Region + ".matbin");

 // Lake
 IntegerVector Lake_cellNumber_int = load_vecbin(path_HydroParam + "Lake_cellNumber_int" + "_" + name_Region + ".vecbin");
 NumericVector Lake_area_km2 = load_vecbin(path_HydroParam + "Lake_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector Lake_capacity_m3 = load_vecbin(path_HydroParam + "Lake_capacity_m3" + "_" + name_Region + ".vecbin");
 NumericVector Lake_albedo_1 = load_vecbin(path_HydroParam + "Lake_albedo_1" + "_" + name_Region + ".vecbin");

 // Riverlak
 IntegerVector Riverlak_cellNumber_int = load_vecbin(path_HydroParam + "Riverlak_cellNumber_int" + "_" + name_Region + ".vecbin");
 NumericVector Riverlak_area_km2 = load_vecbin(path_HydroParam + "Riverlak_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector Riverlak_capacity_m3 = load_vecbin(path_HydroParam + "Riverlak_capacity_m3" + "_" + name_Region + ".vecbin");
 NumericVector Riverlak_albedo_1 = load_vecbin(path_HydroParam + "Riverlak_albedo_1" + "_" + name_Region + ".vecbin");

 // Reservoi
 IntegerVector Reservoi_cellNumber_int = load_vecbin(path_HydroParam + "Reservoi_cellNumber_int" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_area_km2 = load_vecbin(path_HydroParam + "Reservoi_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_capacity_m3 = load_vecbin(path_HydroParam + "Reservoi_capacity_m3" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_albedo_1 = load_vecbin(path_HydroParam + "Reservoi_albedo_1" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_meanDemand_m3 = load_vecbin(path_HydroParam + "Reservoi_meanDemand_m3" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_meanInflow_m3 = load_vecbin(path_HydroParam + "Reservoi_meanInflow_m3" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_releaseCoefficient_1 = load_vecbin(path_HydroParam + "Reservoi_releaseCoefficient_1" + "_" + name_Region + ".vecbin");
 NumericVector Reservoi_dayIrrigate_int = load_vecbin(path_HydroParam + "Reservoi_dayIrrigate_int" + "_" + name_Region + ".vecbin");

 // // demand
 NumericMatrix Reservoi_demand_m3_READ = load_matbin(path_WateruseInput + "Reservoi_demand_m3" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix Reservoi_demand_m3(idx_Month.size(), Reservoi_demand_m3_READ.ncol());
 idx_Month = idx_Month - 1;
 for (int i = 0; i < idx_Month.size(); ++i) {
   Reservoi_demand_m3.row(i) = Reservoi_demand_m3_READ.row(idx_Month[i]);
 }

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
   SNOW_ice_mm = load_vecbin(path_InitialState + "SNOW_ice_mm" + "_" + name_Region + ".vecbin");
   LAND_interceptWater_mm = load_vecbin(path_InitialState + "LAND_interceptWater_mm" + "_" + name_Region + ".vecbin");
   SOIL_water_mm = load_vecbin(path_InitialState + "SOIL_water_mm" + "_" + name_Region + ".vecbin");
   GROUND_water_mm = load_vecbin(path_InitialState + "GROUND_water_mm" + "_" + name_Region + ".vecbin");
   RIVER_water_m3 = load_vecbin(path_InitialState + "RIVER_water_m3" + "_" + name_Region + ".vecbin");
   Lake_water_m3 = load_vecbin(path_InitialState + "Lake_water_m3" + "_" + name_Region + ".vecbin");
   Riverlak_water_m3 = load_vecbin(path_InitialState + "Riverlak_water_m3" + "_" + name_Region + ".vecbin");
 }

 // Parameters
 NumericVector param_ATMOS_thr_Ts = load_vecbin(path_HydroParam + "param_ATMOS_thr_Ts" + "_" + name_Region + ".vecbin");
 NumericVector param_SNOW_fac_f = load_vecbin(path_HydroParam + "param_SNOW_fac_f" + "_" + name_Region + ".vecbin");
 NumericVector param_SNOW_fac_Tmelt = load_vecbin(path_HydroParam + "param_SNOW_fac_Tmelt" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_prt_alpha = load_vecbin(path_HydroParam + "param_EVATRANS_prt_alpha" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_vic_gamma = load_vecbin(path_HydroParam + "param_EVATRANS_vic_gamma" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_sup_k = load_vecbin(path_HydroParam + "param_EVATRANS_sup_k" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_sup_gamma = load_vecbin(path_HydroParam + "param_EVATRANS_sup_gamma" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_wat_petmax = load_vecbin(path_HydroParam + "param_EVATRANS_wat_petmax" + "_" + name_Region + ".vecbin");
 NumericVector param_INFILT_hbv_beta = load_vecbin(path_HydroParam + "param_INFILT_hbv_beta" + "_" + name_Region + ".vecbin");
 LogicalVector param_PERCOLA_wat_01 = load_vecbin(path_HydroParam + "param_PERCOLA_wat_01" + "_" + name_Region + ".vecbin");
 NumericVector param_PERCOLA_wat_k = load_vecbin(path_HydroParam + "param_PERCOLA_wat_k" + "_" + name_Region + ".vecbin");
 NumericVector param_PERCOLA_wat_thresh = load_vecbin(path_HydroParam + "param_PERCOLA_wat_thresh" + "_" + name_Region + ".vecbin");
 NumericVector param_BASEFLOW_sur_k = load_vecbin(path_HydroParam + "param_BASEFLOW_sur_k" + "_" + name_Region + ".vecbin");
 NumericVector param_Lake_acp_storeFactor = load_vecbin(path_HydroParam + "param_Lake_acp_storeFactor" + "_" + name_Region + ".vecbin");
 NumericVector param_Lake_acp_gamma = load_vecbin(path_HydroParam + "param_Lake_acp_gamma" + "_" + name_Region + ".vecbin");
 NumericVector param_Riverlak_lin_storeFactor = load_vecbin(path_HydroParam + "param_Riverlak_lin_storeFactor" + "_" + name_Region + ".vecbin");

 // Call the main WaterGAP3_N function (this would need to be implemented in C++ as well)





 NumericMatrix CELL_discharge_m3 = WaterGAP3_U(name_Region,
                                             mark_Time,
                                             n_Time,
                                             n_Spat,
                                             ATMOS_precipitation_mm,
                                             ATMOS_temperature_Cel,
                                             ATMOS_solarRadiat_MJ,
                                             ATMOS_solarRadiatClearSky_MJ,
                                             WITHDRAW_surface_m3,
                                             WITHDRAW_ground_m3,
                                             Reservoi_demand_m3,
                                             CELL_cellNumberAround_int,
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
                                             Reservoi_cellNumber_int,
                                             Reservoi_area_km2,
                                             Reservoi_capacity_m3,
                                             Reservoi_albedo_1,
                                             Reservoi_meanDemand_m3,
                                             Reservoi_meanInflow_m3,
                                             Reservoi_releaseCoefficient_1,
                                             Reservoi_dayIrrigate_int,
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


//' @rdname run_model
//' @export
// [[Rcpp::export]]
NumericMatrix run_WaterGAP3_N(std::string name_Region, std::string mark_Time, int n_Time, int n_Spat,
                             std::string path_MeteoInput, std::string path_HydroParam,
                             std::string path_InitialState = "UNKNOW",
                             std::string path_Boundary = "UNKNOW",
                             std::string path_VariExport = "NonExport",
                             std::string path_FinalState = "NonExport") {

 // Meteo Forcing Matrix
 NumericMatrix ATMOS_precipitation_mm = load_matbin(path_MeteoInput + "ATMOS_precipitation_mm" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix ATMOS_temperature_Cel = load_matbin(path_MeteoInput + "ATMOS_temperature_Cel" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix ATMOS_solarRadiat_MJ = load_matbin(path_MeteoInput + "ATMOS_solarRadiat_MJ" + "_" + name_Region + "_" + mark_Time + ".matbin");
 NumericMatrix ATMOS_solarRadiatClearSky_MJ = load_matbin(path_MeteoInput + "ATMOS_solarRadiatClearSky_MJ" + "_" + name_Region + "_" + mark_Time + ".matbin");

 // Calculate LAND_interceptCapacity_mm
 NumericMatrix LAND_leafAreaIndex_1 = load_matbin(path_MeteoInput + "LAND_leafAreaIndex_1" + "_" + name_Region + "_" + mark_Time + ".matbin");
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
   Upstream_cellNumber_int = load_vecbin(path_Boundary + "Upstream_cellNumber_int" + "_" + name_Region + ".vecbin");
   Upstream_streamflow_m3 = load_matbin(path_Boundary + "Upstream_streamflow_m3" + "_" + name_Region + ".matbin");
 }
 // Hydro Parameter Vector
 // LAND
 NumericVector LAND_area_km2 = load_vecbin(path_HydroParam + "LAND_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector LAND_albedo_1 = load_vecbin(path_HydroParam + "LAND_albedo_1" + "_" + name_Region + ".vecbin");
 NumericVector LAND_snowAlbedo_1 = load_vecbin(path_HydroParam + "LAND_snowAlbedo_1" + "_" + name_Region + ".vecbin");
 NumericVector LAND_builtRatio_1 = load_vecbin(path_HydroParam + "LAND_builtRatio_1" + "_" + name_Region + ".vecbin");

 // SOIL
 NumericVector SOIL_capacity_mm = load_vecbin(path_HydroParam + "SOIL_capacity_mm" + "_" + name_Region + ".vecbin");
 NumericVector SOIL_potentialPercola_mm = load_vecbin(path_HydroParam + "SOIL_potentialPercola_mm" + "_" + name_Region + ".vecbin");

 // RIVER
 NumericVector RIVER_length_km = load_vecbin(path_HydroParam + "RIVER_length_km" + "_" + name_Region + ".vecbin");
 NumericVector RIVER_velocity_km = load_vecbin(path_HydroParam + "RIVER_velocity_km" + "_" + name_Region + ".vecbin");

 // CELL
 NumericVector CELL_elevation_m = load_vecbin(path_HydroParam + "CELL_elevation_m" + "_" + name_Region + ".vecbin");
 List CELL_cellNumberStep_int = read_int_vector_list(path_HydroParam + "CELL_cellNumberStep_int" + "_" + name_Region + ".bin");
 List CELL_inflowCellNumberStep_int = read_int_matrix_list(path_HydroParam + "CELL_inflowCellNumberStep_int" + "_" + name_Region + ".bin");

 // Lake
 IntegerVector Lake_cellNumber_int = load_vecbin(path_HydroParam + "Lake_cellNumber_int" + "_" + name_Region + ".vecbin");
 NumericVector Lake_area_km2 = load_vecbin(path_HydroParam + "Lake_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector Lake_capacity_m3 = load_vecbin(path_HydroParam + "Lake_capacity_m3" + "_" + name_Region + ".vecbin");
 NumericVector Lake_albedo_1 = load_vecbin(path_HydroParam + "Lake_albedo_1" + "_" + name_Region + ".vecbin");

 // Riverlak
 IntegerVector Riverlak_cellNumber_int = load_vecbin(path_HydroParam + "Riverlak_cellNumber_int" + "_" + name_Region + ".vecbin");
 NumericVector Riverlak_area_km2 = load_vecbin(path_HydroParam + "Riverlak_area_km2" + "_" + name_Region + ".vecbin");
 NumericVector Riverlak_capacity_m3 = load_vecbin(path_HydroParam + "Riverlak_capacity_m3" + "_" + name_Region + ".vecbin");
 NumericVector Riverlak_albedo_1 = load_vecbin(path_HydroParam + "Riverlak_albedo_1" + "_" + name_Region + ".vecbin");

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
   SNOW_ice_mm = load_vecbin(path_InitialState + "SNOW_ice_mm" + "_" + name_Region + ".vecbin");
   LAND_interceptWater_mm = load_vecbin(path_InitialState + "LAND_interceptWater_mm" + "_" + name_Region + ".vecbin");
   SOIL_water_mm = load_vecbin(path_InitialState + "SOIL_water_mm" + "_" + name_Region + ".vecbin");
   GROUND_water_mm = load_vecbin(path_InitialState + "GROUND_water_mm" + "_" + name_Region + ".vecbin");
   RIVER_water_m3 = load_vecbin(path_InitialState + "RIVER_water_m3" + "_" + name_Region + ".vecbin");
   Lake_water_m3 = load_vecbin(path_InitialState + "Lake_water_m3" + "_" + name_Region + ".vecbin");
   Riverlak_water_m3 = load_vecbin(path_InitialState + "Riverlak_water_m3" + "_" + name_Region + ".vecbin");
 }

 // Parameters
 NumericVector param_ATMOS_thr_Ts = load_vecbin(path_HydroParam + "param_ATMOS_thr_Ts" + "_" + name_Region + ".vecbin");
 NumericVector param_SNOW_fac_f = load_vecbin(path_HydroParam + "param_SNOW_fac_f" + "_" + name_Region + ".vecbin");
 NumericVector param_SNOW_fac_Tmelt = load_vecbin(path_HydroParam + "param_SNOW_fac_Tmelt" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_prt_alpha = load_vecbin(path_HydroParam + "param_EVATRANS_prt_alpha" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_vic_gamma = load_vecbin(path_HydroParam + "param_EVATRANS_vic_gamma" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_sup_k = load_vecbin(path_HydroParam + "param_EVATRANS_sup_k" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_sup_gamma = load_vecbin(path_HydroParam + "param_EVATRANS_sup_gamma" + "_" + name_Region + ".vecbin");
 NumericVector param_EVATRANS_wat_petmax = load_vecbin(path_HydroParam + "param_EVATRANS_wat_petmax" + "_" + name_Region + ".vecbin");
 NumericVector param_INFILT_hbv_beta = load_vecbin(path_HydroParam + "param_INFILT_hbv_beta" + "_" + name_Region + ".vecbin");
 LogicalVector param_PERCOLA_wat_01 = load_vecbin(path_HydroParam + "param_PERCOLA_wat_01" + "_" + name_Region + ".UNF1");
 NumericVector param_PERCOLA_wat_k = load_vecbin(path_HydroParam + "param_PERCOLA_wat_k" + "_" + name_Region + ".vecbin");
 NumericVector param_PERCOLA_wat_thresh = load_vecbin(path_HydroParam + "param_PERCOLA_wat_thresh" + "_" + name_Region + ".vecbin");
 NumericVector param_BASEFLOW_sur_k = load_vecbin(path_HydroParam + "param_BASEFLOW_sur_k" + "_" + name_Region + ".vecbin");
 NumericVector param_Lake_acp_storeFactor = load_vecbin(path_HydroParam + "param_Lake_acp_storeFactor" + "_" + name_Region + ".vecbin");
 NumericVector param_Lake_acp_gamma = load_vecbin(path_HydroParam + "param_Lake_acp_gamma" + "_" + name_Region + ".vecbin");
 NumericVector param_Riverlak_lin_storeFactor = load_vecbin(path_HydroParam + "param_Riverlak_lin_storeFactor" + "_" + name_Region + ".vecbin");

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
