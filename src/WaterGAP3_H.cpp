#include "WaterGAP3.h"

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
   NumericVector SNOW_ice_mm,
   NumericVector SOIL_water_mm,
   NumericVector SOIL_capacity_mm,
   NumericVector SOIL_potentialPercola_mm,
   NumericVector GROUND_water_mm,
   NumericVector GROUND_capacity_mm,
   NumericVector RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   NumericVector CELL_landArea_km2,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   NumericVector param_ATMOS_thr_Ts,
   NumericVector param_SNOW_fac_f,
   NumericVector param_SNOW_fac_Tmelt,
   NumericVector param_EVATRANS_tur_k,
   NumericVector param_EVATRANS_ubc_gamma,
   NumericVector param_INFILT_hbv_beta,
   NumericVector param_PERCOLA_arn_k,
   NumericVector param_PERCOLA_arn_thresh,
   NumericVector param_BASEFLOW_grf_gamma,
   bool if_allVariExport = false
)
{

 NumericVector ATMOS_rainFall_mm, ATMOS_snowFall_mm,
               SNOW_melt_mm,
               LAND_water_mm, LAND_runoff_mm,
               SOIL_evatrans_mm, SOIL_infilt_mm, SOIL_percola_mm,
               GROUND_basefloW_mm,
               CELL_outflow_m3;
 NumericMatrix RIVER_outflow_m3(n_time, n_spat);

 NumericMatrix OUT_snow(n_time, n_spat), OUT_evatrans(n_time, n_spat),
 OUT_landrunoff, OUT_groundbaseflow,
 OUT_soilwater(n_time, n_spat), OUT_groundwater(n_time, n_spat),
 OUT_snowice(n_time, n_spat), OUT_snowmelt(n_time, n_spat),
 OUT_cellOutflow(n_time, n_spat),
 OUT_riverwater(n_time, n_spat);



 for (int i = 0; i < n_time; i++) {


   // Vertical
   CELL_outflow_m3 = module_land_WaterGAP3(
     ATMOS_precipitation_mm(i, _),
     ATMOS_temperature_Cel(i, _),
     ATMOS_potentialEvatrans_mm(i, _),
     ATMOS_snowFall_mm,
     SNOW_ice_mm,
     LAND_runoff_mm,
     SOIL_water_mm,
     SOIL_capacity_mm,
     SOIL_potentialPercola_mm,
     SOIL_evatrans_mm,
     GROUND_water_mm,
     GROUND_capacity_mm,
     GROUND_basefloW_mm,
     CELL_landArea_km2,
     param_ATMOS_thr_Ts,
     param_EVATRANS_ubc_gamma,
     param_SNOW_fac_f,
     param_SNOW_fac_Tmelt,
     param_INFILT_hbv_beta,
     param_PERCOLA_arn_thresh,
     param_PERCOLA_arn_k,
     param_BASEFLOW_grf_gamma);

   // Horizontal
   RIVER_outflow_m3(i, _) = confluen_WaterGAP3(
     CELL_outflow_m3,
     RIVER_water_m3,
     RIVER_length_km,
     RIVER_velocity_km,
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
     _["snowFall_mm"] = OUT_snow,
     _["evatrans_mm"] = OUT_evatrans,
     _["soilwater_mm"] = OUT_soilwater,
     _["groundwater_mm"] = OUT_groundwater,
     _["snowice_mm"] = OUT_snowice,
     _["runoff_mm"] = OUT_landrunoff,
     _["basefloW_mm"] = OUT_groundbaseflow,
     _["riverwater_m3"] = OUT_riverwater,
     _["streamflow_m3"] = RIVER_outflow_m3
   );
 } else{
   return List::create(
     _["streamflow_m3"] = RIVER_outflow_m3
   );
 }


}
