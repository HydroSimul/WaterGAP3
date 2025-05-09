#include "module.h"
// [[Rcpp::interfaces(r, cpp)]]


//' model
//' @name model
//' @inheritParams HydroGallery::all_vari
//' @param name_Region A short identifier of continents: "eu", "af, "as", "au", "na", "sa".
//' @param mark_Time A short identifier of time.
//' @param path_VariExport Directory path where model output variables (e.g., soil moisture, runoff) will be saved. If "NonExport", no variables will be exported.
//' @param path_FinalState Directory path where final state variables will be saved. If "NonExport", final states will not be saved.
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
NumericMatrix WaterGAP3_U(
    std::string name_Region,
    std::string mark_Time,
    int n_time,
    int n_spat,
    NumericMatrix ATMOS_precipitation_mm,
    NumericMatrix ATMOS_temperature_Cel,
    NumericMatrix ATMOS_solarRadiat_MJ,
    NumericMatrix ATMOS_solarRadiatClearSky_MJ,
    NumericMatrix WITHDRAW_surface_m3,
    NumericMatrix WITHDRAW_ground_m3,
    NumericMatrix Reservoi_demand_m3,
    IntegerMatrix CELL_cellNumberAround_int,
    IntegerVector Upstream_cellNumber_int,
    NumericMatrix Upstream_streamflow_m3,
    NumericVector SNOW_ice_mm,
    NumericVector LAND_area_km2,
    NumericVector LAND_albedo_1,
    NumericVector LAND_snowAlbedo_1,
    NumericVector LAND_builtRatio_1,
    NumericVector LAND_interceptWater_mm,
    NumericMatrix LAND_interceptCapacity_mm,
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialPercola_mm,
    NumericVector GROUND_water_mm,
    NumericVector RIVER_water_m3,
    NumericVector RIVER_length_km,
    NumericVector RIVER_velocity_km,
    NumericVector CELL_elevation_m,
    List CELL_cellNumberStep_int,
    List CELL_inflowCellNumberStep_int,
    IntegerVector Lake_cellNumber_int,
    NumericVector Lake_water_m3,
    NumericVector Lake_area_km2,
    NumericVector Lake_capacity_m3,
    NumericVector Lake_albedo_1,
    IntegerVector Riverlak_cellNumber_int,
    NumericVector Riverlak_water_m3,
    NumericVector Riverlak_area_km2,
    NumericVector Riverlak_capacity_m3,
    NumericVector Riverlak_albedo_1,
    IntegerVector Reservoi_cellNumber_int,
    NumericVector Reservoi_area_km2,
    NumericVector Reservoi_capacity_m3,
    NumericVector Reservoi_albedo_1,
    NumericVector Reservoi_meanInflow_m3,
    NumericVector Reservoi_meanDemand_m3,
    NumericVector Reservoi_releaseCoefficient_1,
    NumericVector Reservoi_dayIrrigate_int,
    NumericVector param_ATMOS_thr_Ts,
    NumericVector param_SNOW_fac_f,
    NumericVector param_SNOW_fac_Tmelt,
    NumericVector param_EVATRANS_prt_alpha,
    NumericVector param_EVATRANS_vic_gamma,
    NumericVector param_EVATRANS_sup_k,
    NumericVector param_EVATRANS_sup_gamma,
    NumericVector param_EVATRANS_wat_petmax,
    NumericVector param_INFILT_hbv_beta,
    LogicalVector param_PERCOLA_wat_01,
    NumericVector param_PERCOLA_wat_k,
    NumericVector param_PERCOLA_wat_thresh,
    NumericVector param_BASEFLOW_sur_k,
    NumericVector param_Lake_acp_storeFactor,
    NumericVector param_Lake_acp_gamma,
    NumericVector param_Riverlak_lin_storeFactor,
    std::string path_FinalState = "NonExport",
    std::string path_VariExport = "NonExport"
)
{

 NumericVector ATMOS_rainFall_mm, ATMOS_snowFall_mm, ATMOS_potentialEvatrans_mm,
               SNOW_melt_mm,
               LAND_water_mm, LAND_runoff_mm,
               SOIL_evatrans_mm, SOIL_infilt_mm, SOIL_percola_mm,
               GROUND_basefloW_mm,
               CELL_innerOutflow_m3;
 NumericVector CELL_withdrawSurface_m3(n_spat), CELL_withdrawGround_mm(n_spat);
 NumericMatrix RIVER_outflow_m3(n_time, n_spat);

 if (Upstream_cellNumber_int(0) != 0) {
   int n_Upstream = Upstream_streamflow_m3.ncol();
   for (int i = 0; i < n_Upstream; i++) {
     RIVER_outflow_m3(_, Upstream_cellNumber_int[i] - 1) = Upstream_streamflow_m3(_, i);
   }
 }

 NumericMatrix OUT_snow(n_time, n_spat), OUT_evatrans(n_time, n_spat), OUT_potentialevatrans(n_time, n_spat),
 OUT_landrunoff(n_time, n_spat), OUT_groundbaseflow(n_time, n_spat),
 OUT_soilwater(n_time, n_spat), OUT_groundwater(n_time, n_spat),
 OUT_snowice(n_time, n_spat), OUT_snowmelt(n_time, n_spat),
 OUT_cellOutflow(n_time, n_spat),
 OUT_riverwater(n_time, n_spat);

 NumericVector Lake_evatrans_mm, Riverlak_evatrans_mm, Reservoi_evatrans_mm;
 int n_Lake = Lake_cellNumber_int.size(), n_Riverlak = Riverlak_cellNumber_int.size(), n_Reservoi = Reservoi_cellNumber_int.size();
 NumericMatrix OUT_lakeWater(n_time, n_Lake),
 OUT_lakeEvalake(n_time, n_Lake), OUT_riverlakEvalake(n_time, n_Riverlak), OUT_reservoiEvalake(n_time, n_Reservoi);




 for (int i = 0; i < n_time; i++) {
  LogicalVector Reservoi_isIrrigate_01 = (Reservoi_dayIrrigate_int == i);


   // Vertical
   CELL_innerOutflow_m3 = module_land_WaterGAP3(
     ATMOS_precipitation_mm(i, _),
     ATMOS_temperature_Cel(i, _),
     ATMOS_solarRadiat_MJ(i, _),
     ATMOS_solarRadiatClearSky_MJ(i, _),
     ATMOS_potentialEvatrans_mm,
     ATMOS_snowFall_mm,
     SNOW_ice_mm,
     LAND_area_km2,
     LAND_albedo_1,
     LAND_snowAlbedo_1,
     LAND_builtRatio_1,
     LAND_interceptWater_mm,
     LAND_interceptCapacity_mm(i, _),
     LAND_runoff_mm,
     SOIL_water_mm,
     SOIL_capacity_mm,
     SOIL_potentialPercola_mm,
     SOIL_evatrans_mm,
     GROUND_water_mm,
     GROUND_basefloW_mm,
     CELL_elevation_m,
     param_ATMOS_thr_Ts,
     param_EVATRANS_prt_alpha,
     param_EVATRANS_sup_k,
     param_EVATRANS_sup_gamma,
     param_EVATRANS_wat_petmax,
     param_SNOW_fac_f,
     param_SNOW_fac_Tmelt,
     param_INFILT_hbv_beta,
     param_PERCOLA_wat_01,
     param_PERCOLA_wat_thresh,
     param_PERCOLA_wat_k,
     param_BASEFLOW_sur_k);



   // Local Lake
   NumericVector Lake_verticalInflow_m3 = module_waterbody_WaterGAP3(
     subset_get(ATMOS_precipitation_mm(i, _), Lake_cellNumber_int),
     subset_get(ATMOS_temperature_Cel(i, _), Lake_cellNumber_int),
     subset_get(ATMOS_solarRadiat_MJ(i, _), Lake_cellNumber_int),
     subset_get(ATMOS_solarRadiatClearSky_MJ(i, _), Lake_cellNumber_int),
     Lake_water_m3,
     Lake_area_km2,
     Lake_capacity_m3,
     Lake_evatrans_mm,
     Lake_albedo_1,
     subset_get(CELL_elevation_m, Lake_cellNumber_int),
     subset_get(param_EVATRANS_prt_alpha, Lake_cellNumber_int),
     subset_get(param_EVATRANS_vic_gamma, Lake_cellNumber_int));


   NumericVector Lake_Outflow_m3 = module_lake_WaterGAP3(
     Lake_water_m3,
     Lake_capacity_m3,
     Lake_verticalInflow_m3,
     subset_get(CELL_innerOutflow_m3, Lake_cellNumber_int),
     param_Lake_acp_storeFactor,
     param_Lake_acp_gamma);


   subset_put(CELL_innerOutflow_m3, Lake_cellNumber_int, Lake_Outflow_m3);


   // Riverlak
   NumericVector Riverlak_water_m3 = subset_get(RIVER_water_m3, Riverlak_cellNumber_int);
   NumericVector Riverlak_verticalInflow_m3 = module_waterbody_WaterGAP3(
     subset_get(ATMOS_precipitation_mm(i, _), Riverlak_cellNumber_int),
     subset_get(ATMOS_temperature_Cel(i, _), Riverlak_cellNumber_int),
     subset_get(ATMOS_solarRadiat_MJ(i, _), Riverlak_cellNumber_int),
     subset_get(ATMOS_solarRadiatClearSky_MJ(i, _), Riverlak_cellNumber_int),
     Riverlak_water_m3,
     Riverlak_area_km2,
     Riverlak_capacity_m3,
     Riverlak_evatrans_mm,
     Riverlak_albedo_1,
     subset_get(CELL_elevation_m, Riverlak_cellNumber_int),
     subset_get(param_EVATRANS_prt_alpha, Riverlak_cellNumber_int),
     subset_get(param_EVATRANS_vic_gamma, Riverlak_cellNumber_int));
   Riverlak_water_m3 = pmin(Riverlak_water_m3 + Riverlak_verticalInflow_m3, Riverlak_capacity_m3);
   subset_put(RIVER_water_m3, Riverlak_cellNumber_int, Riverlak_water_m3);
   // when vertical inflow more than lake capacity, the overflow will added into Cell-innerOutflow
   subset_put(CELL_innerOutflow_m3, Riverlak_cellNumber_int, pmax(Riverlak_water_m3 + Riverlak_verticalInflow_m3 - Riverlak_capacity_m3,0.));
   // Reservoi
   NumericVector Reservoi_water_m3 = subset_get(RIVER_water_m3, Reservoi_cellNumber_int);
   NumericVector Reservoi_verticalInflow_m3 = module_waterbody_WaterGAP3(
     subset_get(ATMOS_precipitation_mm(i, _), Reservoi_cellNumber_int),
     subset_get(ATMOS_temperature_Cel(i, _), Reservoi_cellNumber_int),
     subset_get(ATMOS_solarRadiat_MJ(i, _), Reservoi_cellNumber_int),
     subset_get(ATMOS_solarRadiatClearSky_MJ(i, _), Reservoi_cellNumber_int),
     Reservoi_water_m3,
     Reservoi_area_km2,
     Reservoi_capacity_m3,
     Reservoi_evatrans_mm,
     Reservoi_albedo_1,
     subset_get(CELL_elevation_m, Reservoi_cellNumber_int),
     subset_get(param_EVATRANS_prt_alpha, Reservoi_cellNumber_int),
     subset_get(param_EVATRANS_vic_gamma, Reservoi_cellNumber_int));
   Reservoi_water_m3 = pmin(Reservoi_water_m3 + Reservoi_verticalInflow_m3, Reservoi_capacity_m3);
   subset_put(RIVER_water_m3, Reservoi_cellNumber_int, Reservoi_water_m3);
   // when vertical inflow more than lake capacity, the overflow will added into Cell-innerOutflow
   subset_put(CELL_innerOutflow_m3, Reservoi_cellNumber_int, pmax(Reservoi_water_m3 + Reservoi_verticalInflow_m3 - Reservoi_capacity_m3,0.));



   // Horizontal
   // Wateruse withdrawal
   CELL_withdrawGround_mm += WITHDRAW_ground_m3(i, _) / LAND_area_km2 / 1000;
   withdraw_SingleCell(CELL_withdrawGround_mm,
                          GROUND_water_mm);

   CELL_withdrawSurface_m3 += WITHDRAW_surface_m3(i, _);
   withdrawSurface_Around(CELL_withdrawSurface_m3,
                          RIVER_water_m3,
                          Lake_water_m3,
                          CELL_cellNumberAround_int);

   // confluen
   if (Upstream_cellNumber_int(0) != 0) {
     subset_put(CELL_innerOutflow_m3, Upstream_cellNumber_int, Upstream_streamflow_m3(i,_));
   }
   NumericVector RIVER_outflow_m3_TEMP = confluen_WaterGAP3(
     RIVER_water_m3,
     RIVER_length_km,
     RIVER_velocity_km,
     CELL_innerOutflow_m3,
     CELL_cellNumberStep_int,
     CELL_inflowCellNumberStep_int,
     Riverlak_cellNumber_int,
     Riverlak_capacity_m3,
     Reservoi_cellNumber_int,
     Reservoi_demand_m3(i, _),
     Reservoi_capacity_m3,
     Reservoi_meanInflow_m3,
     Reservoi_meanDemand_m3,
     Reservoi_releaseCoefficient_1,
     Reservoi_isIrrigate_01,
     param_Riverlak_lin_storeFactor
   );

   RIVER_outflow_m3(i, _) = RIVER_outflow_m3_TEMP;


   if (path_VariExport != "NonExport") {
     OUT_potentialevatrans(i, _) = ATMOS_potentialEvatrans_mm;
     OUT_snow(i, _) = ATMOS_snowFall_mm;
     OUT_evatrans(i, _) = SOIL_evatrans_mm;
     OUT_landrunoff(i, _) = LAND_runoff_mm;
     OUT_soilwater(i, _) = SOIL_water_mm;
     OUT_groundwater(i, _) = GROUND_water_mm;
     OUT_groundbaseflow(i, _) = GROUND_basefloW_mm;
     OUT_snowice(i, _) = SNOW_ice_mm;
     OUT_riverwater(i, _) = RIVER_water_m3;
     OUT_lakeWater(i, _) = Lake_water_m3;
     OUT_lakeEvalake(i, _) = Lake_evatrans_mm;
     OUT_riverlakEvalake(i, _) = Riverlak_evatrans_mm;
     OUT_reservoiEvalake(i, _) = Reservoi_evatrans_mm;
   }

 }

 if (path_FinalState != "NonExport") {
   save_vecbin(SNOW_ice_mm,            path_FinalState + "SNOW_ice_mm_"            + name_Region + ".vecbin");
   save_vecbin(LAND_interceptWater_mm, path_FinalState + "LAND_interceptWater_mm_" + name_Region + ".vecbin");
   save_vecbin(SOIL_water_mm,          path_FinalState + "SOIL_water_mm_"          + name_Region + ".vecbin");
   save_vecbin(GROUND_water_mm,        path_FinalState + "GROUND_water_mm_"        + name_Region + ".vecbin");
   save_vecbin(RIVER_water_m3,         path_FinalState + "RIVER_water_m3_"         + name_Region + ".vecbin");
   save_vecbin(Lake_water_m3,          path_FinalState + "Lake_water_m3_"          + name_Region + ".vecbin");
 }

 if (path_VariExport != "NonExport") {
   save_matbin(RIVER_outflow_m3,     path_VariExport + "CELL_discharge_m3_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_potentialevatrans,path_VariExport + "ATMOS_potentialEvatrans_mm_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_snow,             path_VariExport + "ATMOS_snowFall_mm_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_evatrans,         path_VariExport + "SOIL_evatrans_mm_"     + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_landrunoff,       path_VariExport + "LAND_runoff_mm_"       + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_soilwater,        path_VariExport + "SOIL_water_mm_"        + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_groundwater,      path_VariExport + "GROUND_water_mm_"      + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_groundbaseflow,   path_VariExport + "GROUND_baseflow_mm_"   + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_snowice,          path_VariExport + "SNOW_ice_mm_"          + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_riverwater,       path_VariExport + "RIVER_water_m3_"       + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_lakeWater,        path_VariExport + "Lake_water_m3_"        + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_lakeEvalake,      path_VariExport + "Lake_evatrans_mm_"     + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_riverlakEvalake,  path_VariExport + "Riverlak_evatrans_mm_" + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_reservoiEvalake,  path_VariExport + "Seservoi_evatrans_mm_" + name_Region + "_" + mark_Time + ".matbin");
 }

 return RIVER_outflow_m3;


}






//' @rdname model
//' @export
// [[Rcpp::export]]
NumericMatrix WaterGAP3_N(
   std::string name_Region,
   std::string mark_Time,
   int n_time,
   int n_spat,
   NumericMatrix ATMOS_precipitation_mm,
   NumericMatrix ATMOS_temperature_Cel,
   NumericMatrix ATMOS_solarRadiat_MJ,
   NumericMatrix ATMOS_solarRadiatClearSky_MJ,
   IntegerVector Upstream_cellNumber_int,
   NumericMatrix Upstream_streamflow_m3,
   NumericVector SNOW_ice_mm,
   NumericVector LAND_area_km2,
   NumericVector LAND_albedo_1,
   NumericVector LAND_snowAlbedo_1,
   NumericVector LAND_builtRatio_1,
   NumericVector LAND_interceptWater_mm,
   NumericMatrix LAND_interceptCapacity_mm,
   NumericVector SOIL_water_mm,
   NumericVector SOIL_capacity_mm,
   NumericVector SOIL_potentialPercola_mm,
   NumericVector GROUND_water_mm,
   NumericVector RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   NumericVector CELL_elevation_m,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   IntegerVector Lake_cellNumber_int,
   NumericVector Lake_water_m3,
   NumericVector Lake_area_km2,
   NumericVector Lake_capacity_m3,
   NumericVector Lake_albedo_1,
   IntegerVector Riverlak_cellNumber_int,
   NumericVector Riverlak_area_km2,
   NumericVector Riverlak_capacity_m3,
   NumericVector Riverlak_albedo_1,
   NumericVector param_ATMOS_thr_Ts,
   NumericVector param_SNOW_fac_f,
   NumericVector param_SNOW_fac_Tmelt,
   NumericVector param_EVATRANS_prt_alpha,
   NumericVector param_EVATRANS_vic_gamma,
   NumericVector param_EVATRANS_sup_k,
   NumericVector param_EVATRANS_sup_gamma,
   NumericVector param_EVATRANS_wat_petmax,
   NumericVector param_INFILT_hbv_beta,
   LogicalVector param_PERCOLA_wat_01,
   NumericVector param_PERCOLA_wat_k,
   NumericVector param_PERCOLA_wat_thresh,
   NumericVector param_BASEFLOW_sur_k,
   NumericVector param_Lake_acp_storeFactor,
   NumericVector param_Lake_acp_gamma,
   NumericVector param_Riverlak_lin_storeFactor,
   std::string path_FinalState = "NonExport",
   std::string path_VariExport = "NonExport"
)
{

 NumericVector ATMOS_rainFall_mm, ATMOS_snowFall_mm, ATMOS_potentialEvatrans_mm,
 SNOW_melt_mm,
 LAND_water_mm, LAND_runoff_mm,
 SOIL_evatrans_mm, SOIL_infilt_mm, SOIL_percola_mm,
 GROUND_basefloW_mm,
 CELL_innerOutflow_m3;
 NumericMatrix RIVER_outflow_m3(n_time, n_spat);

 if (Upstream_cellNumber_int(0) != 0) {
   int n_Upstream = Upstream_streamflow_m3.ncol();
   for (int i = 0; i < n_Upstream; i++) {
     RIVER_outflow_m3(_, Upstream_cellNumber_int[i] - 1) = Upstream_streamflow_m3(_, i);
   }
 }

 NumericMatrix OUT_snow(n_time, n_spat), OUT_evatrans(n_time, n_spat), OUT_potentialevatrans(n_time, n_spat),
 OUT_landrunoff(n_time, n_spat), OUT_groundbaseflow(n_time, n_spat),
 OUT_soilwater(n_time, n_spat), OUT_groundwater(n_time, n_spat),
 OUT_snowice(n_time, n_spat), OUT_snowmelt(n_time, n_spat),
 OUT_cellOutflow(n_time, n_spat),
 OUT_riverwater(n_time, n_spat);

 NumericVector Lake_evatrans_mm, Riverlak_evatrans_mm;
 int n_Lake = Lake_cellNumber_int.size(), n_Riverlak = Riverlak_cellNumber_int.size();
 NumericMatrix OUT_lakeWater(n_time, n_Lake), OUT_riverlakWater(n_time, n_Riverlak),
 OUT_lakeEvalake(n_time, n_Lake), OUT_riverlakEvalake(n_time, n_Riverlak);




 for (int i = 0; i < n_time; i++) {


   // Vertical
   CELL_innerOutflow_m3 = module_land_WaterGAP3(
     ATMOS_precipitation_mm(i, _),
     ATMOS_temperature_Cel(i, _),
     ATMOS_solarRadiat_MJ(i, _),
     ATMOS_solarRadiatClearSky_MJ(i, _),
     ATMOS_potentialEvatrans_mm,
     ATMOS_snowFall_mm,
     SNOW_ice_mm,
     LAND_area_km2,
     LAND_albedo_1,
     LAND_snowAlbedo_1,
     LAND_builtRatio_1,
     LAND_interceptWater_mm,
     LAND_interceptCapacity_mm(i, _),
     LAND_runoff_mm,
     SOIL_water_mm,
     SOIL_capacity_mm,
     SOIL_potentialPercola_mm,
     SOIL_evatrans_mm,
     GROUND_water_mm,
     GROUND_basefloW_mm,
     CELL_elevation_m,
     param_ATMOS_thr_Ts,
     param_EVATRANS_prt_alpha,
     param_EVATRANS_sup_k,
     param_EVATRANS_sup_gamma,
     param_EVATRANS_wat_petmax,
     param_SNOW_fac_f,
     param_SNOW_fac_Tmelt,
     param_INFILT_hbv_beta,
     param_PERCOLA_wat_01,
     param_PERCOLA_wat_thresh,
     param_PERCOLA_wat_k,
     param_BASEFLOW_sur_k);



   // Local Lake
   NumericVector Lake_verticalInflow_m3 = module_waterbody_WaterGAP3(
     subset_get(ATMOS_precipitation_mm(i, _), Lake_cellNumber_int),
     subset_get(ATMOS_temperature_Cel(i, _), Lake_cellNumber_int),
     subset_get(ATMOS_solarRadiat_MJ(i, _), Lake_cellNumber_int),
     subset_get(ATMOS_solarRadiatClearSky_MJ(i, _), Lake_cellNumber_int),
     Lake_water_m3,
     Lake_area_km2,
     Lake_capacity_m3,
     Lake_evatrans_mm,
     Lake_albedo_1,
     subset_get(CELL_elevation_m, Lake_cellNumber_int),
     subset_get(param_EVATRANS_prt_alpha, Lake_cellNumber_int),
     subset_get(param_EVATRANS_vic_gamma, Lake_cellNumber_int));

   NumericVector Lake_Outflow_m3 = module_lake_WaterGAP3(
     Lake_water_m3,
     Lake_capacity_m3,
     Lake_verticalInflow_m3,
     subset_get(CELL_innerOutflow_m3, Lake_cellNumber_int),
     param_Lake_acp_storeFactor,
     param_Lake_acp_gamma);


   subset_put(CELL_innerOutflow_m3, Lake_cellNumber_int, Lake_Outflow_m3);


   // Riverlak
   NumericVector Riverlak_water_m3 = subset_get(RIVER_water_m3, Riverlak_cellNumber_int);
   NumericVector Riverlak_verticalInflow_m3 = module_waterbody_WaterGAP3(
     subset_get(ATMOS_precipitation_mm(i, _), Riverlak_cellNumber_int),
     subset_get(ATMOS_temperature_Cel(i, _), Riverlak_cellNumber_int),
     subset_get(ATMOS_solarRadiat_MJ(i, _), Riverlak_cellNumber_int),
     subset_get(ATMOS_solarRadiatClearSky_MJ(i, _), Riverlak_cellNumber_int),
     Riverlak_water_m3,
     Riverlak_area_km2,
     Riverlak_capacity_m3,
     Riverlak_evatrans_mm,
     Riverlak_albedo_1,
     subset_get(CELL_elevation_m, Riverlak_cellNumber_int),
     subset_get(param_EVATRANS_prt_alpha, Riverlak_cellNumber_int),
     subset_get(param_EVATRANS_vic_gamma, Riverlak_cellNumber_int));
   Riverlak_water_m3 = pmin(Riverlak_water_m3 + Riverlak_verticalInflow_m3, Riverlak_capacity_m3);
   subset_put(RIVER_water_m3, Riverlak_cellNumber_int, Riverlak_water_m3);
   // when vertical inflow more than lake capacity, the overflow will added into Cell-innerOutflow
   subset_put(CELL_innerOutflow_m3, Riverlak_cellNumber_int, pmax(Riverlak_water_m3 + Riverlak_verticalInflow_m3 - Riverlak_capacity_m3,0.));


   // confluen
   // Horizontal
   if (Upstream_cellNumber_int(0) != 0) {
     subset_put(CELL_innerOutflow_m3, Upstream_cellNumber_int, Upstream_streamflow_m3(i,_));
   }
   NumericVector RIVER_outflow_m3_TEMP = confluen_WaterGAP3_N(
     RIVER_water_m3,
     RIVER_length_km,
     RIVER_velocity_km,
     CELL_innerOutflow_m3, // just for the cell in Step 0
     CELL_cellNumberStep_int,
     CELL_inflowCellNumberStep_int,
     Riverlak_cellNumber_int,
     Riverlak_capacity_m3,
     param_Riverlak_lin_storeFactor
   );

   RIVER_outflow_m3(i, _) = RIVER_outflow_m3_TEMP;


   if (path_VariExport != "NonExport") {
     OUT_potentialevatrans(i, _) = ATMOS_potentialEvatrans_mm;
     OUT_snow(i, _) = ATMOS_snowFall_mm;
     OUT_evatrans(i, _) = SOIL_evatrans_mm;
     OUT_landrunoff(i, _) = LAND_runoff_mm;
     OUT_soilwater(i, _) = SOIL_water_mm;
     OUT_groundwater(i, _) = GROUND_water_mm;
     OUT_groundbaseflow(i, _) = GROUND_basefloW_mm;
     OUT_snowice(i, _) = SNOW_ice_mm;
     OUT_riverwater(i, _) = RIVER_water_m3;
     OUT_lakeWater(i, _) = Lake_water_m3;
     OUT_lakeEvalake(i, _) = Lake_evatrans_mm;
     // OUT_riverlakWater(i, _) = Riverlak_water_m3;
     OUT_riverlakEvalake(i, _) = Riverlak_evatrans_mm;
   }

 }

 if (path_FinalState != "NonExport") {
   save_vecbin(SNOW_ice_mm,            path_FinalState + "SNOW_ice_mm_"            + name_Region + ".vecbin");
   save_vecbin(LAND_interceptWater_mm, path_FinalState + "LAND_interceptWater_mm_" + name_Region + ".vecbin");
   save_vecbin(SOIL_water_mm,          path_FinalState + "SOIL_water_mm_"          + name_Region + ".vecbin");
   save_vecbin(GROUND_water_mm,        path_FinalState + "GROUND_water_mm_"        + name_Region + ".vecbin");
   save_vecbin(RIVER_water_m3,         path_FinalState + "RIVER_water_m3_"         + name_Region + ".vecbin");
   save_vecbin(Lake_water_m3,          path_FinalState + "Lake_water_m3_"          + name_Region + ".vecbin");
   // save_vecbin(Riverlak_water_m3,      path_FinalState + "Riverlak_water_m3_"      + name_Region + ".vecbin");
 }

 if (path_VariExport != "NonExport") {
   save_matbin(RIVER_outflow_m3,     path_VariExport + "CELL_discharge_m3_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_potentialevatrans,path_VariExport + "ATMOS_potentialEvatrans_mm_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_snow,             path_VariExport + "ATMOS_snowFall_mm_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_evatrans,         path_VariExport + "SOIL_evatrans_mm_"     + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_landrunoff,       path_VariExport + "LAND_runoff_mm_"       + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_soilwater,        path_VariExport + "SOIL_water_mm_"        + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_groundwater,      path_VariExport + "GROUND_water_mm_"      + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_groundbaseflow,   path_VariExport + "GROUND_baseflow_mm_"   + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_snowice,          path_VariExport + "SNOW_ice_mm_"          + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_riverwater,       path_VariExport + "RIVER_water_m3_"       + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_lakeWater,        path_VariExport + "Lake_water_m3_"        + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_lakeEvalake,      path_VariExport + "Lake_evatrans_mm_"     + name_Region + "_" + mark_Time + ".matbin");
   // save_matbin(OUT_riverlakWater,    path_VariExport + "Riverlak_water_m3_"    + name_Region + "_" + mark_Time + ".matbin");
   save_matbin(OUT_riverlakEvalake,  path_VariExport + "Riverlak_evatrans_mm_" + name_Region + "_" + mark_Time + ".matbin");
 }

 return RIVER_outflow_m3;


}
