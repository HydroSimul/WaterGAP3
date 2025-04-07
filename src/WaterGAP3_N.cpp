#include "WaterGAP3_Module.h"
// [[Rcpp::interfaces(r, cpp)]]

//' WaterGAP3
//' @name WaterGAP3
//' @inheritParams HydroGallery::all_vari
//' @param name_Region A short identifier of continents: "eu", "af, "as", "au", "na", "sa".
//' @param path_VariExport Directory path (as a string) where model output variables (e.g., soil moisture, runoff) will be saved.
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
NumericMatrix WaterGAP3_N(
    std::string name_Region,
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
    NumericVector Riverlak_water_m3,
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

 NumericVector ATMOS_rainFall_mm, ATMOS_snowFall_mm,
               SNOW_melt_mm,
               LAND_water_mm, LAND_runoff_mm,
               SOIL_evatrans_mm, SOIL_infilt_mm, SOIL_percola_mm,
               GROUND_basefloW_mm,
               CELL_verticalflow_m3;
 NumericMatrix RIVER_outflow_m3(n_time, n_spat);

 if (Upstream_cellNumber_int(0) != 0) {
   int n_Upstream = Upstream_streamflow_m3.ncol();
   for (int i = 0; i < n_Upstream; i++) {
     RIVER_outflow_m3(_, Upstream_cellNumber_int[i] - 1) = Upstream_streamflow_m3(_, i);
   }
 }

 NumericMatrix OUT_snow(n_time, n_spat), OUT_evatrans(n_time, n_spat),
 OUT_landrunoff(n_time, n_spat), OUT_groundbaseflow(n_time, n_spat),
 OUT_soilwater(n_time, n_spat), OUT_groundwater(n_time, n_spat),
 OUT_snowice(n_time, n_spat), OUT_snowmelt(n_time, n_spat),
 OUT_cellOutflow(n_time, n_spat),
 OUT_riverwater(n_time, n_spat);

 NumericVector Lake_evatrans_mm, Riverlak_evatrans_mm;
 int n_Lake = Lake_cellNumber_int.size(), n_Riverlake = Riverlak_cellNumber_int.size();
 NumericMatrix OUT_lakeWater(n_time, n_Lake), OUT_riverlakWater(n_time, n_Riverlake),
 OUT_lakeEvalake(n_time, n_Lake), OUT_riverlakEvalake(n_time, n_Riverlake);




 for (int i = 0; i < n_time; i++) {


   // Vertical
   CELL_verticalflow_m3 = module_land_WaterGAP3(
     ATMOS_precipitation_mm(i, _),
     ATMOS_temperature_Cel(i, _),
     ATMOS_solarRadiat_MJ(i, _),
     ATMOS_solarRadiatClearSky_MJ(i, _),
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
     subset_get(CELL_verticalflow_m3, Lake_cellNumber_int),
     param_Lake_acp_storeFactor,
     param_Lake_acp_gamma);


   subset_put(CELL_verticalflow_m3, Lake_cellNumber_int, Lake_Outflow_m3);


   // Riverlake
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

   Riverlak_water_m3 = pmax(Riverlak_water_m3 + subset_get(CELL_verticalflow_m3, Riverlak_cellNumber_int) + Riverlak_verticalInflow_m3, 0.);


   // confluen
   RIVER_water_m3 += CELL_verticalflow_m3;

   // Horizontal
   if (Upstream_cellNumber_int(0) != 0) {
     subset_put(CELL_verticalflow_m3, Upstream_cellNumber_int, Upstream_streamflow_m3(i,_));
   }
   RIVER_outflow_m3(i, _) = confluen_WaterGAP3_L(
     RIVER_water_m3,
     RIVER_length_km,
     RIVER_velocity_km,
     CELL_verticalflow_m3, // just for the cell in Step 0
     CELL_cellNumberStep_int,
     CELL_inflowCellNumberStep_int,
     Riverlak_cellNumber_int,
     Riverlak_water_m3,
     Riverlak_capacity_m3,
     param_Riverlak_lin_storeFactor
   );

   if (path_VariExport != "NonExport") {
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
     OUT_riverlakWater(i, _) = Riverlak_water_m3;
     OUT_riverlakEvalake(i, _) = Riverlak_evatrans_mm;
   }

 }

 if (path_FinalState != "NonExport") {
   write_unf(SNOW_ice_mm,            path_FinalState + "SNOW_ice_mm_"            + name_Region + ".UNF0");
   write_unf(LAND_interceptWater_mm, path_FinalState + "LAND_interceptWater_mm_" + name_Region + ".UNF0");
   write_unf(SOIL_water_mm,          path_FinalState + "SOIL_water_mm_"          + name_Region + ".UNF0");
   write_unf(GROUND_water_mm,        path_FinalState + "GROUND_water_mm_"        + name_Region + ".UNF0");
   write_unf(RIVER_water_m3,         path_FinalState + "RIVER_water_m3_"         + name_Region + ".UNF0");
   write_unf(Lake_water_m3,          path_FinalState + "Lake_water_m3_"          + name_Region + ".UNF0");
   write_unf(Riverlak_water_m3,      path_FinalState + "Riverlak_water_m3_"      + name_Region + ".UNF0");
 }

 if (path_VariExport != "NonExport") {
   save_wgmat(RIVER_outflow_m3,     path_VariExport + "CELL_discharge_m3_"    + name_Region + ".wgmat");
   save_wgmat(OUT_snow,             path_VariExport + "ATMOS_snowFall_mm_"    + name_Region + ".wgmat");
   save_wgmat(OUT_evatrans,         path_VariExport + "SOIL_evatrans_mm_"     + name_Region + ".wgmat");
   save_wgmat(OUT_landrunoff,       path_VariExport + "LAND_runoff_mm_"       + name_Region + ".wgmat");
   save_wgmat(OUT_soilwater,        path_VariExport + "SOIL_water_mm_"        + name_Region + ".wgmat");
   save_wgmat(OUT_groundwater,      path_VariExport + "GROUND_water_mm_"      + name_Region + ".wgmat");
   save_wgmat(OUT_groundbaseflow,   path_VariExport + "GROUND_baseflow_mm_"   + name_Region + ".wgmat");
   save_wgmat(OUT_snowice,          path_VariExport + "SNOW_ice_mm_"          + name_Region + ".wgmat");
   save_wgmat(OUT_riverwater,       path_VariExport + "RIVER_water_m3_"       + name_Region + ".wgmat");
   save_wgmat(OUT_lakeWater,        path_VariExport + "Lake_water_m3_"        + name_Region + ".wgmat");
   save_wgmat(OUT_lakeEvalake,      path_VariExport + "Lake_evatrans_mm_"     + name_Region + ".wgmat");
   save_wgmat(OUT_riverlakWater,    path_VariExport + "Riverlak_water_m3_"    + name_Region + ".wgmat");
   save_wgmat(OUT_riverlakEvalake,  path_VariExport + "Riverlak_evatrans_mm_" + name_Region + ".wgmat");
 }

 return RIVER_outflow_m3;


}
