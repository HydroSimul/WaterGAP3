#include "WaterGAP3_Process.h"
// [[Rcpp::interfaces(r, cpp)]]

//' WaterGAP3
//' @name WaterGAP3_H
//' @inheritParams all_vari
//' @inheritParams process
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
List WaterGAP3_H(
   int n_time,
   int n_spat,
   NumericMatrix ATMOS_precipitation_mm,
   NumericMatrix ATMOS_temperature_Cel,
   NumericMatrix ATMOS_potentialEvatrans_mm,
   IntegerVector Upstream_cellNumber_int,
   NumericMatrix Upstream_streamflow_m3,
   NumericVector SNOW_ice_mm,
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
   NumericVector CELL_landArea_km2,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   NumericVector param_ATMOS_thr_Ts,
   NumericVector param_SNOW_fac_f,
   NumericVector param_SNOW_fac_Tmelt,
   NumericVector param_EVATRANS_sup_k,
   NumericVector param_EVATRANS_sup_gamma,
   NumericVector param_EVATRANS_wat_petmax,
   NumericVector param_INFILT_hbv_beta,
   LogicalVector param_PERCOLA_wat_01,
   NumericVector param_PERCOLA_wat_k,
   NumericVector param_PERCOLA_wat_thresh,
   NumericVector param_BASEFLOW_sur_k,
   bool if_allVariExport = false
)
{

 NumericVector ATMOS_rainFall_mm, ATMOS_snowFall_mm,
               SNOW_melt_mm,
               LAND_water_mm, LAND_runoff_mm,
               SOIL_evatrans_mm, SOIL_infilt_mm, SOIL_percola_mm,
               GROUND_basefloW_mm,
               CELL_verticalflow_m3;
 NumericMatrix RIVER_outflow_m3(n_time, n_spat);

 NumericMatrix OUT_snow(n_time, n_spat), OUT_evatrans(n_time, n_spat),
 OUT_landrunoff(n_time, n_spat), OUT_groundbaseflow(n_time, n_spat),
 OUT_soilwater(n_time, n_spat), OUT_groundwater(n_time, n_spat),
 OUT_snowice(n_time, n_spat), OUT_snowmelt(n_time, n_spat),
 OUT_cellOutflow(n_time, n_spat),
 OUT_riverwater(n_time, n_spat);



 for (int i = 0; i < n_time; i++) {


   // Vertical
   CELL_verticalflow_m3 = module_land_WaterGAP3(
     ATMOS_precipitation_mm(i, _),
     ATMOS_temperature_Cel(i, _),
     ATMOS_potentialEvatrans_mm(i, _),
     ATMOS_snowFall_mm,
     SNOW_ice_mm,
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
     CELL_landArea_km2,
     param_ATMOS_thr_Ts,
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

   RIVER_water_m3 += CELL_verticalflow_m3;

   // Horizontal
   if (Upstream_cellNumber_int(0) != 0) {
     subset_put(CELL_verticalflow_m3, Upstream_cellNumber_int, Upstream_streamflow_m3(i,_));
   }

   RIVER_outflow_m3(i, _) = confluen_WaterGAP3(
     RIVER_water_m3,
     RIVER_length_km,
     RIVER_velocity_km,
     CELL_verticalflow_m3, // just for the cell in Step 0
     CELL_cellNumberStep_int,
     CELL_inflowCellNumberStep_int
   );

   if (if_allVariExport) {
     OUT_snow(i, _) = ATMOS_snowFall_mm;
     OUT_evatrans(i, _) = SOIL_evatrans_mm;
     OUT_landrunoff(i, _) = LAND_runoff_mm;
     OUT_soilwater(i, _) = SOIL_water_mm;
     OUT_groundwater(i, _) = GROUND_water_mm;
     OUT_groundbaseflow(i, _) = GROUND_basefloW_mm;
     OUT_snowice(i, _) = SNOW_ice_mm;
     OUT_riverwater(i, _) = RIVER_water_m3;
   }

 }

 if (if_allVariExport) {
   return List::create(
     _["streamflow_m3"] = RIVER_outflow_m3,
     // _["snowFall_mm"] = OUT_snow,
     // _["evatrans_mm"] = OUT_evatrans,
     // _["soilwater_mm"] = OUT_soilwater,
     // _["groundwater_mm"] = OUT_groundwater,
     // _["snowice_mm"] = OUT_snowice,
     // _["runoff_mm"] = OUT_landrunoff,
     // _["basefloW_mm"] = OUT_groundbaseflow,
     _["riverwater_m3"] = OUT_riverwater
   );
 } else{
   return List::create(
     _["streamflow_m3"] = RIVER_outflow_m3
   );
 }


}
