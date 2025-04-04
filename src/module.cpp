#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

NumericVector module_land_WaterGAP3(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector ATMOS_solarRadiatClearSky_MJ,
    NumericVector& ATMOS_snowFall_mm,
    NumericVector& SNOW_ice_mm,
    NumericVector LAND_albedo_1,
    NumericVector LAND_snowAlbedo_1,
    NumericVector LAND_builtRatio_1,
    NumericVector& LAND_interceptWater_mm,
    NumericVector LAND_interceptCapacity_mm,
    NumericVector& LAND_runoff_mm,
    NumericVector& SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialPercola_mm,
    NumericVector& SOIL_evatrans_mm,
    NumericVector& GROUND_water_mm,
    NumericVector& GROUND_basefloW_mm,
    NumericVector CELL_landArea_km2,
    NumericVector CELL_elevation_m,
    NumericVector param_ATMOS_thr_Ts,
    NumericVector param_EVATRANS_prt_alpha,
    NumericVector param_EVATRANS_sup_k,
    NumericVector param_EVATRANS_sup_gamma,
    NumericVector param_EVATRANS_wat_petmax,
    NumericVector param_SNOW_fac_f,
    NumericVector param_SNOW_fac_Tmelt,
    NumericVector param_INFILT_hbv_beta,
    LogicalVector param_PERCOLA_wat_01,
    NumericVector param_PERCOLA_wat_thresh,
    NumericVector param_PERCOLA_wat_k,
    NumericVector param_BASEFLOW_sur_k) {

  NumericVector ATMOS_rainFall_mm, SNOW_melt_mm, LAND_water_mm, SOIL_INFILT_mm, SOIL_percola_mm;
  // // snow fall
  ATMOS_snowFall_mm = atmosSnow_ThresholdT(ATMOS_precipitation_mm, ATMOS_temperature_Cel, param_ATMOS_thr_Ts);
  ATMOS_rainFall_mm = ATMOS_precipitation_mm - ATMOS_snowFall_mm;

  // // Interception
  NumericVector LAND_intercp = intercep_Full(ATMOS_rainFall_mm, LAND_interceptWater_mm, LAND_interceptCapacity_mm);
  LAND_interceptWater_mm += LAND_intercp;
  ATMOS_rainFall_mm += -LAND_intercp;

  // // Potential Evapo
  LAND_albedo_1 = ifelse(SNOW_ice_mm > 3, LAND_snowAlbedo_1, LAND_albedo_1);
  NumericVector ATMOS_netRadiat_MJ = meteo_nettoRadiat_WaterGAP3(
    ATMOS_temperature_Cel,
    ATMOS_solarRadiat_MJ,
    ATMOS_solarRadiatClearSky_MJ,
    LAND_albedo_1);
  NumericVector ATMOS_potentialEvatrans_mm = evatransPotential_PriestleyTaylor(
    ATMOS_temperature_Cel,
    ATMOS_netRadiat_MJ,
    CELL_elevation_m,
    param_EVATRANS_prt_alpha);
  // // Evapo Interception
  NumericVector LAND_intercepEvapo_mm = evatransActual_SupplyPow(ATMOS_potentialEvatrans_mm,
                                                                 LAND_interceptWater_mm, LAND_interceptCapacity_mm,
                                                                 param_EVATRANS_sup_k, param_EVATRANS_sup_gamma);
  LAND_interceptWater_mm += -LAND_intercepEvapo_mm;
  ATMOS_potentialEvatrans_mm += -LAND_intercepEvapo_mm;

  LAND_water_mm = ATMOS_rainFall_mm;
  // // Built up
  NumericVector LAND_runoffBuiltup_mm = LAND_water_mm * (LAND_builtRatio_1 * .5);
  LAND_water_mm += -LAND_runoffBuiltup_mm;


  // // Evapo soil
  SOIL_evatrans_mm = evatransActual_WaterGAP3(ATMOS_potentialEvatrans_mm, SOIL_water_mm, SOIL_capacity_mm, param_EVATRANS_wat_petmax);
  SOIL_water_mm += -SOIL_evatrans_mm;


  SOIL_evatrans_mm += LAND_intercepEvapo_mm;
  // // Snow melt
  SNOW_melt_mm = snowMelt_Factor(SNOW_ice_mm, ATMOS_temperature_Cel, param_SNOW_fac_f, param_SNOW_fac_Tmelt);
  LAND_water_mm += SNOW_melt_mm;
  SNOW_ice_mm += -SNOW_melt_mm;
  SNOW_ice_mm += ATMOS_snowFall_mm;

  // // soil infiltration
  SOIL_INFILT_mm = infilt_HBV(LAND_water_mm, SOIL_water_mm, SOIL_capacity_mm, param_INFILT_hbv_beta);
  SOIL_water_mm += SOIL_INFILT_mm;
  LAND_runoff_mm = LAND_water_mm - SOIL_INFILT_mm;

  // // percolation
  SOIL_percola_mm = percola_WaterGAP3(LAND_runoff_mm, SOIL_potentialPercola_mm, param_PERCOLA_wat_01, param_PERCOLA_wat_thresh, param_PERCOLA_wat_k);
  GROUND_water_mm += SOIL_percola_mm;
  LAND_runoff_mm += LAND_runoffBuiltup_mm - SOIL_percola_mm;

  // // ground water
  GROUND_basefloW_mm = baseflow_SupplyRatio(GROUND_water_mm, param_BASEFLOW_sur_k);
  GROUND_water_mm += - GROUND_basefloW_mm;



  return (LAND_runoff_mm + GROUND_basefloW_mm) * CELL_landArea_km2 * 1000;
}

NumericVector module_waterbody_WaterGAP3(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector ATMOS_solarRadiatClearSky_MJ,
    NumericVector Lake_water_m3,
    NumericVector Lake_area_km2,
    NumericVector Lake_capacity_m3,
    NumericVector& Lake_evatrans_mm,
    NumericVector Lake_albedo_1,
    NumericVector CELL_elevation_m,
    NumericVector param_EVATRANS_prt_alpha,
    NumericVector param_EVATRANS_vic_gamma) {
  // // Potential Evapo
  NumericVector ATMOS_netRadiat_MJ = meteo_nettoRadiat_WaterGAP3(
    ATMOS_temperature_Cel,
    ATMOS_solarRadiat_MJ,
    ATMOS_solarRadiatClearSky_MJ,
    Lake_albedo_1);
  NumericVector ATMOS_potentialEvatrans_mm = evatransPotential_PriestleyTaylor(
    ATMOS_temperature_Cel,
    ATMOS_netRadiat_MJ,
    CELL_elevation_m,
    param_EVATRANS_prt_alpha);

  // Step 1: Calculate Lake evaporation in mm
  Lake_evatrans_mm = evatransActual_VIC(
    ATMOS_potentialEvatrans_mm,
    Lake_water_m3,
    Lake_capacity_m3,
    param_EVATRANS_vic_gamma
  );

  // Step 2: Calculate vertical inflow in mm and convert to m3
  return (ATMOS_precipitation_mm - Lake_evatrans_mm) * Lake_area_km2 * 1000; // can be negativ

}



NumericVector module_lake_WaterGAP3(
    NumericVector& Lake_water_m3,
    NumericVector Lake_capacity_m3,
    NumericVector Lake_verticalInflow_m3,
    NumericVector Lake_inflow_m3,
    NumericVector param_Lake_acp_storeFactor,
    NumericVector param_Lake_acp_gamma) {




  // Step 3: Update Lake water storage with inflows
  Lake_water_m3 = pmax(Lake_water_m3 + Lake_inflow_m3 + Lake_verticalInflow_m3, 0);
  NumericVector Lake_overflow_m3 = pmax(Lake_water_m3 -  Lake_capacity_m3, 0);
  Lake_water_m3 = pmin(Lake_water_m3, Lake_capacity_m3);

  // Step 4: Calculate Lake outflow in m3
  NumericVector Lake_Outflow_m3 = lakeout_SupplyPow(
    Lake_water_m3,
    Lake_capacity_m3,
    param_Lake_acp_storeFactor,
    param_Lake_acp_gamma
  );

  // Step 5: Update Lake water storage with outflows and enforce capacity constraints
  Lake_water_m3 += -Lake_Outflow_m3;
  // Lake_water_m3 = pmin(Lake_water_m3, Lake_capacity_m3);
  // Return updated Lake water storage
  return Lake_Outflow_m3 + Lake_overflow_m3;
}


















