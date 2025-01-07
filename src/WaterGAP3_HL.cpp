#include "WaterGAP3.h"

//' WaterGAP3
//' @name WaterGAP3_HL
//' @inheritParams all_vari
//' @inheritParams process
//' @param param_Evalake_vic_gamma same as param_EVATRANS_vic_gamma
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
List WaterGAP3_HL(
   int n_time,
   int n_spat,
   NumericMatrix ATMOS_precipitation_mm,
   NumericMatrix ATMOS_temperature_Cel,
   NumericMatrix ATMOS_potentialEvatrans_mm,
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
   IntegerVector Lake_cellNumber_int,
   NumericVector Lake_water_m3,
   NumericVector Lake_area_km2,
   NumericVector Lake_capacity_m3,
   IntegerVector Riverlak_cellNumber_int,
   NumericVector Riverlak_water_m3,
   NumericVector Riverlak_area_km2,
   NumericVector Riverlak_capacity_m3,
   NumericVector param_ATMOS_thr_Ts,
   NumericVector param_SNOW_fac_f,
   NumericVector param_SNOW_fac_Tmelt,
   NumericVector param_EVATRANS_sup_k,
   NumericVector param_EVATRANS_sup_gamma,
   NumericVector param_EVATRANS_sur_k,
   NumericVector param_INFILT_hbv_beta,
   NumericVector param_PERCOLA_arn_k,
   NumericVector param_PERCOLA_arn_thresh,
   NumericVector param_BASEFLOW_sur_k,
   NumericVector param_Evalake_vic_gamma,
   NumericVector param_Lake_acp_storeFactor,
   NumericVector param_Lake_acp_gamma,
   NumericVector param_Riverlak_lin_storeFactor,
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
   CELL_outflow_m3 = module_land_WaterGAP3(
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
     param_EVATRANS_sur_k,
     param_SNOW_fac_f,
     param_SNOW_fac_Tmelt,
     param_INFILT_hbv_beta,
     param_PERCOLA_arn_thresh,
     param_PERCOLA_arn_k,
     param_BASEFLOW_sur_k);


   // AET Waterbody
   Lake_evatrans_mm = evatransActual_VIC(
     subset_get(ATMOS_potentialEvatrans_mm, Lake_cellNumber_int),
     Lake_water_m3,
     Lake_capacity_m3,
     subset_get(param_Evalake_vic_gamma, Lake_cellNumber_int)
   );
   Riverlak_evatrans_mm = evatransActual_VIC(
     subset_get(ATMOS_potentialEvatrans_mm, Riverlak_cellNumber_int),
     Riverlak_water_m3,
     Riverlak_capacity_m3,
     subset_get(param_Evalake_vic_gamma, Riverlak_cellNumber_int)
   );


   // // Net vertical Inflow
   NumericVector Lake_verticalInflow_m3, Riverlak_verticalInflow_m3, reservoir_verticalInflow_m3;

   NumericVector Lake_verticalInflow_mm = pmax((subset_get(ATMOS_precipitation_mm, Lake_cellNumber_int) - Lake_evatrans_mm), 0.0);
   Lake_verticalInflow_m3 = Lake_verticalInflow_mm * Lake_area_km2 * 1000;

   NumericVector Riverlak_verticalInflow_mm = pmax((subset_get(ATMOS_precipitation_mm, Riverlak_cellNumber_int) - Riverlak_evatrans_mm), 0.0);
   Riverlak_verticalInflow_m3 = Riverlak_verticalInflow_mm * Riverlak_area_km2 * 1000;

   // Local Lake runoff
   NumericVector Lake_Outflow_m3 = lake_AcceptPow(
     Lake_water_m3,
     (subset_get(CELL_outflow_m3, Lake_cellNumber_int) + Lake_verticalInflow_m3),
     Lake_capacity_m3,
     subset_get(param_Lake_acp_storeFactor, Lake_cellNumber_int),
     subset_get(param_Lake_acp_gamma, Lake_cellNumber_int)
   );

   // Sum Inflow to water net
   subset_put(CELL_outflow_m3, Lake_cellNumber_int, Lake_Outflow_m3);
   subset_put(CELL_outflow_m3, Riverlak_cellNumber_int, subset_get(CELL_outflow_m3, Riverlak_cellNumber_int) + Riverlak_verticalInflow_m3);

   // Horizontal
   RIVER_outflow_m3(i, _) = confluen_WaterGAP3_L(
     CELL_outflow_m3,
     RIVER_water_m3,
     RIVER_length_km,
     RIVER_velocity_km,
     CELL_cellNumberStep_int,
     CELL_inflowCellNumberStep_int,
     Riverlak_cellNumber_int,
     Riverlak_water_m3,
     Riverlak_capacity_m3,
     param_Riverlak_lin_storeFactor
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
     OUT_lakeWater(i, _) = Lake_water_m3;
     OUT_lakeEvalake(i, _) = Lake_evatrans_mm;
     OUT_riverlakWater(i, _) = Riverlak_water_m3;
     OUT_riverlakEvalake(i, _) = Riverlak_evatrans_mm;
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
     _["streamflow_m3"] = RIVER_outflow_m3,
     _["lakewater_m3"] = OUT_lakeWater,
     _["lakeEva_mm"] = OUT_lakeEvalake,
     _["riverlakwater_m3"] = OUT_riverlakWater,
     _["riverlakEva_mm"] = OUT_riverlakEvalake
   );
 } else{
   return List::create(
     _["streamflow_m3"] = RIVER_outflow_m3
   );
 }


}
