#include "WaterGAP3_Module.h"
// [[Rcpp::interfaces(r, cpp)]]

//' WaterGAP3
//' @name WaterGAP3
//' @inheritParams HydroGallery::all_vari
//' @param param_Evalake_vic_gamma same as param_EVATRANS_vic_gamma
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
List WaterGAP3_N(
   int n_time,
   int n_spat,
   NumericMatrix ATMOS_precipitation_mm,
   NumericMatrix ATMOS_temperature_Cel,
   NumericMatrix ATMOS_solarRadiat_MJ,
   NumericMatrix ATMOS_solarRadiatClearSky_MJ,
   IntegerVector Upstream_cellNumber_int,
   NumericMatrix Upstream_streamflow_m3,
   NumericVector SNOW_ice_mm,
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
   NumericVector CELL_landArea_km2,
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
   NumericVector param_Lake_Eva_vic_gamma,
   NumericVector param_Lake_acp_storeFactor,
   NumericVector param_Lake_acp_gamma,
   NumericVector param_Riverlak_Eva_vic_gamma,
   NumericVector param_Riverlak_lin_storeFactor,
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
     CELL_landArea_km2,
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
     CELL_elevation_m,
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
     CELL_elevation_m,
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
     _["streamflow_m3"] = RIVER_outflow_m3,
     _["snowFall_mm"] = OUT_snow,
     _["evatrans_mm"] = OUT_evatrans,
     _["soilwater_mm"] = OUT_soilwater,
     _["groundwater_mm"] = OUT_groundwater,
     _["snowice_mm"] = OUT_snowice,
     _["runoff_mm"] = OUT_landrunoff,
     _["baseflow_mm"] = OUT_groundbaseflow,
     _["riverwater_m3"] = OUT_riverwater,
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
