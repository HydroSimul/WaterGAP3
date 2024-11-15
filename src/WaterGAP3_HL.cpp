#include "WaterGAP3.h"

//' WaterGAP3
//' @name WaterGAP3
//' @inheritParams all_vari
//' @inheritParams process
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
List WaterGAP3_HLR(
   int n_time,
   int n_spat,
   NumericMatrix atmos_precipitation_mm,
   NumericMatrix atmos_temperature_Cel,
   NumericMatrix atmos_solarRadiat_MJ,
   NumericVector snow_ice_mm,
   NumericVector land_impermeableFrac_1,
   NumericVector soil_water_mm,
   NumericVector soil_capacity_mm,
   NumericVector soil_potentialPercola_mm,
   NumericVector ground_water_mm,
   NumericVector ground_capacity_mm,
   NumericVector basin_landArea_km2,
   IntegerVector lake_cellNumber_int,
   NumericVector lake_water_m3,
   NumericVector lake_area_km2,
   NumericVector lake_capacity_m3,
   NumericVector river_water_m3,
   NumericVector river_length_km,
   NumericVector river_velocity_km,
   IntegerVector riverlake_cellNumber_int,
   NumericVector riverlake_water_m3,
   NumericVector riverlake_area_km2,
   NumericVector riverlake_capacity_m3,
   List basin_cellNumberStep_int,
   List basin_inflowCellNumberStep_int,
   NumericVector param_atmos_thr_Ts,
   NumericVector param_snow_fac_f,
   NumericVector param_snow_fac_Tmelt,
   NumericVector param_evatrans_tur_k,
   NumericVector param_evatrans_ubc_gamma,
   NumericVector param_evatrans_vic_gamma,
   NumericVector param_infilt_hbv_beta,
   NumericVector param_percola_arn_k,
   NumericVector param_percola_arn_thresh,
   NumericVector param_baseflow_grf_gamma,
   NumericVector param_lake_acp_storeFactor,
   NumericVector param_lake_acp_gamma,
   NumericVector param_riverlake_lin_storeFactor
)
{

 NumericVector atmos_potentialEvatrans_mm, atmos_rain_mm, atmos_snow_mm,
 snow_melt_mm, land_water_mm,
 soil_evatrans_mm, soil_infilt_mm, soil_percol_mm,
 lake_evatrans_mm, riverlake_evatrans_mm,
 waterbody_verticalInflow_m3(n_spat, 0.),
 land_outflow_m3;
 NumericMatrix land_runoff_mm(n_time, n_spat), ground_baseflow_mm(n_time, n_spat), river_outflow_m3(n_time, n_spat);
 int n_Lake = lake_cellNumber_int.size(), n_Riverlake = riverlake_cellNumber_int.size();
 NumericMatrix out_snow(n_time, n_spat), out_evatransPot(n_time, n_spat), out_evatrans(n_time, n_spat),
               out_soilwater(n_time, n_spat), out_groundwater(n_time, n_spat),
               out_snowice(n_time, n_spat), out_snowmelt(n_time, n_spat),
               out_land_outflow(n_time, n_spat),
               out_river_water(n_time, n_spat), out_riverlake_water(n_time, n_Riverlake),
               out_lake_evatrans(n_time, n_Lake), out_riverlake_evatrans(n_time, n_Riverlake);


 for (int i = 0; i < n_time; i++) {

   // Land
   // // snow
   atmos_snow_mm = atmosSnow_ThresholdT(atmos_precipitation_mm(i, _), atmos_temperature_Cel(i, _), param_atmos_thr_Ts);
   atmos_rain_mm = atmos_precipitation_mm(i, _) - atmos_snow_mm;

   // // PET
   atmos_potentialEvatrans_mm = evatransPotential_TurcWendling(atmos_temperature_Cel(i, _), atmos_solarRadiat_MJ(i, _), param_evatrans_tur_k);

   // // soil
   soil_evatrans_mm = evatransActual_UBC(atmos_potentialEvatrans_mm, soil_water_mm, soil_capacity_mm, param_evatrans_ubc_gamma);
   soil_water_mm += - soil_evatrans_mm;
   land_water_mm = atmos_rain_mm;

   // // Snow melt
   snow_melt_mm = snowMelt_Factor(snow_ice_mm, atmos_temperature_Cel(i, _), param_snow_fac_f, param_snow_fac_Tmelt);
   land_water_mm += snow_melt_mm;
   snow_ice_mm += -snow_melt_mm;
   snow_ice_mm += atmos_snow_mm;

   // // soil infiltration
   soil_infilt_mm = infilt_HBV(land_water_mm, soil_water_mm, soil_capacity_mm, param_infilt_hbv_beta);
   soil_water_mm += soil_infilt_mm;
   land_runoff_mm(i, _) = land_water_mm - soil_infilt_mm;

   // // soil percolation
   soil_percol_mm = percola_Arno(soil_water_mm, soil_capacity_mm, soil_potentialPercola_mm, param_percola_arn_thresh, param_percola_arn_k);
   ground_water_mm += soil_percol_mm;
   soil_water_mm += - soil_percol_mm;

   // // baseflow
   NumericVector baseflow_temp = ifelse(ground_water_mm < ground_capacity_mm, 0, ground_water_mm - ground_capacity_mm);

   // // ground water
   ground_water_mm = ifelse(ground_water_mm < ground_capacity_mm,ground_water_mm, ground_capacity_mm);
   ground_baseflow_mm(i, _) = baseflow_GR4Jfix(ground_water_mm, ground_capacity_mm, param_baseflow_grf_gamma);
   ground_water_mm += - ground_baseflow_mm(i, _);
   ground_baseflow_mm(i, _) = ground_baseflow_mm(i, _) + baseflow_temp;


   land_outflow_m3 = (land_runoff_mm(i, _) + ground_baseflow_mm(i, _)) * basin_landArea_km2 * 1000;
   // Water body

   // AET Waterbody
   lake_evatrans_mm = evatransActual_VIC(
     subset_get(atmos_potentialEvatrans_mm, lake_cellNumber_int),
     lake_water_m3,
     lake_capacity_m3,
     subset_get(param_evatrans_vic_gamma, lake_cellNumber_int)
   );
   riverlake_evatrans_mm = evatransActual_VIC(
     subset_get(atmos_potentialEvatrans_mm, riverlake_cellNumber_int),
     riverlake_water_m3,
     riverlake_capacity_m3,
     subset_get(param_evatrans_vic_gamma, riverlake_cellNumber_int)
   );

   // // Net vertical Inflow
   NumericVector lake_verticalInflow_m3, riverlake_verticalInflow_m3, reservoir_verticalInflow_m3;
   lake_verticalInflow_m3 = (subset_get(atmos_precipitation_mm, lake_cellNumber_int) - lake_evatrans_mm) * lake_area_km2 * 1000;
   riverlake_verticalInflow_m3 = (subset_get(atmos_precipitation_mm, riverlake_cellNumber_int) - riverlake_evatrans_mm) * riverlake_area_km2 * 1000;

   // Local Lake runoff
   NumericVector lake_Outflow_m3 = lake_AcceptPow(
     lake_water_m3,
     (subset_get(land_outflow_m3, lake_cellNumber_int) + lake_verticalInflow_m3),
     lake_capacity_m3,
     subset_get(param_lake_acp_storeFactor, lake_cellNumber_int),
     subset_get(param_lake_acp_gamma, lake_cellNumber_int)
   );

   // Sum Inflow to water net
   subset_put(land_outflow_m3, lake_cellNumber_int, lake_Outflow_m3);
   subset_put(land_outflow_m3, riverlake_cellNumber_int, subset_get(land_outflow_m3, riverlake_cellNumber_int) + riverlake_verticalInflow_m3);



   river_outflow_m3(i, _) = confluen_WaterGAP3_L(
     land_outflow_m3,
     river_water_m3,
     river_length_km,
     river_velocity_km,
     riverlake_cellNumber_int,
     riverlake_water_m3,
     riverlake_capacity_m3,
     basin_cellNumberStep_int,
     basin_inflowCellNumberStep_int,
     param_lake_acp_storeFactor,
     param_lake_acp_gamma,
     param_riverlake_lin_storeFactor
   );

   // out_snow(i, _) = atmos_snow_mm;
   // out_evatransPot(i, _) = atmos_potentialEvatrans_mm;
   // out_evatrans(i, _) = soil_evatrans_mm;
   // out_soilwater(i, _) = soil_water_mm;
   // out_groundwater(i, _) = ground_water_mm;
   // out_snowice(i, _) = snow_ice_mm;
   // out_snowmelt(i, _) = snow_melt_mm;
   out_land_outflow(i, _) = land_outflow_m3;
   out_river_water(i, _) = river_water_m3;
   out_riverlake_water(i, _) = riverlake_water_m3;
   out_lake_evatrans(i, _) = lake_evatrans_mm;
   out_riverlake_evatrans(i, _) = riverlake_evatrans_mm;

 }

 return List::create(
   // _["precipitation_mm"] = atmos_precipitation_mm,
   // _["temperature_Cel"] = atmos_temperature_Cel,
   // _["solarRadiat_MJ"] = atmos_solarRadiat_MJ,
   // _["snowFall_mm"] = out_snow,
   // _["evatransPot_mm"] = out_evatransPot,
   // _["evatrans_mm"] = out_evatrans,
   // _["soilwater_mm"] = out_soilwater,
   // _["groundwater_mm"] = out_groundwater,
   // _["snowice_mm"] = out_snowice,
   // _["snowmelt_mm"] = out_snowmelt,
   _["runoff_mm"] = land_runoff_mm,
   _["runoff_m3"] = out_land_outflow,
   _["baseflow_mm"] = ground_baseflow_mm,
   _["river_water_m3"] = river_water_m3,
   _["riverlake_water_m3"] = riverlake_water_m3,
   _["lake_evatrans_mm"] = lake_evatrans_mm,
   _["riverlake_evatrans_mm"] = riverlake_evatrans_mm,
   _["streamflow_m3"] = river_outflow_m3
 );
}
