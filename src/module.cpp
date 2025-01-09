#include "WaterGAP3.h"



NumericVector module_land_Sachsen(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector& ATMOS_snowFall_mm,
    NumericVector& SNOW_ice_mm,
    NumericVector& LAND_runoff_mm,
    NumericVector& SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialPercola_mm,
    NumericVector& SOIL_evatrans_mm,
    NumericVector& GROUND_water_mm,
    NumericVector GROUND_capacity_mm,
    NumericVector& GROUND_basefloW_mm,
    NumericVector CELL_landArea_km2,
    NumericVector param_ATMOS_thr_Ts,
    NumericVector param_EVATRANS_ubc_gamma,
    NumericVector param_SNOW_fac_f,
    NumericVector param_SNOW_fac_Tmelt,
    NumericVector param_INFILT_hbv_beta,
    NumericVector param_PERCOLA_arn_thresh,
    NumericVector param_PERCOLA_arn_k,
    NumericVector param_BASEFLOW_grf_gamma) {

  NumericVector ATMOS_rainFall_mm, SNOW_melt_mm, LAND_water_mm, SOIL_INFILT_mm, SOIL_percola_mm;
  // // snow
  ATMOS_snowFall_mm = atmosSnow_ThresholdT(ATMOS_precipitation_mm, ATMOS_temperature_Cel, param_ATMOS_thr_Ts);
  ATMOS_rainFall_mm = ATMOS_precipitation_mm - ATMOS_snowFall_mm;

  // // soil
  SOIL_evatrans_mm = evatransActual_UBC(ATMOS_potentialEvatrans_mm, SOIL_water_mm, SOIL_capacity_mm, param_EVATRANS_ubc_gamma);
  SOIL_water_mm += - SOIL_evatrans_mm;
  LAND_water_mm = ATMOS_rainFall_mm;

  // // Snow melt
  SNOW_melt_mm = snowMelt_Factor(SNOW_ice_mm, ATMOS_temperature_Cel, param_SNOW_fac_f, param_SNOW_fac_Tmelt);
  LAND_water_mm += SNOW_melt_mm;
  SNOW_ice_mm += -SNOW_melt_mm;
  SNOW_ice_mm += ATMOS_snowFall_mm;

  // // soil infiltration
  SOIL_INFILT_mm = infilt_HBV(LAND_water_mm, SOIL_water_mm, SOIL_capacity_mm, param_INFILT_hbv_beta);
  SOIL_water_mm += SOIL_INFILT_mm;
  LAND_runoff_mm = LAND_water_mm - SOIL_INFILT_mm;

  // // soil percolation
  SOIL_percola_mm = percola_Arno(SOIL_water_mm, SOIL_capacity_mm, SOIL_potentialPercola_mm, param_PERCOLA_arn_thresh, param_PERCOLA_arn_k);
  GROUND_water_mm += SOIL_percola_mm;
  SOIL_water_mm += - SOIL_percola_mm;

  // // baseflow
  NumericVector basefloW_temp = ifelse(GROUND_water_mm < GROUND_capacity_mm, 0, GROUND_water_mm - GROUND_capacity_mm);

  // // ground water
  GROUND_water_mm = ifelse(GROUND_water_mm < GROUND_capacity_mm,GROUND_water_mm, GROUND_capacity_mm);
  GROUND_basefloW_mm = baseflow_GR4Jfix(GROUND_water_mm, GROUND_capacity_mm, param_BASEFLOW_grf_gamma);
  GROUND_water_mm += - GROUND_basefloW_mm;
  GROUND_basefloW_mm = GROUND_basefloW_mm + basefloW_temp;



  return (LAND_runoff_mm + GROUND_basefloW_mm) * CELL_landArea_km2 * 1000;
}


NumericVector module_land_WaterGAP3(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector& ATMOS_snowFall_mm,
    NumericVector& SNOW_ice_mm,
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
    NumericVector param_ATMOS_thr_Ts,
    NumericVector param_EVATRANS_sup_k,
    NumericVector param_EVATRANS_sup_gamma,
    NumericVector param_EVATRANS_wat_petmax,
    NumericVector param_SNOW_fac_f,
    NumericVector param_SNOW_fac_Tmelt,
    NumericVector param_INFILT_hbv_beta,
    NumericVector param_PERCOLA_arn_thresh,
    NumericVector param_PERCOLA_arn_k,
    NumericVector param_BASEFLOW_sur_k) {

  NumericVector ATMOS_rainFall_mm, SNOW_melt_mm, LAND_water_mm, SOIL_INFILT_mm, SOIL_percola_mm;
  // // snow
  ATMOS_snowFall_mm = atmosSnow_ThresholdT(ATMOS_precipitation_mm, ATMOS_temperature_Cel, param_ATMOS_thr_Ts);
  ATMOS_rainFall_mm = ATMOS_precipitation_mm - ATMOS_snowFall_mm;

  // // Interception
  NumericVector LAND_intercp = intercep_Full(ATMOS_rainFall_mm, LAND_interceptWater_mm, LAND_interceptCapacity_mm);
  LAND_interceptWater_mm += LAND_intercp;
  ATMOS_rainFall_mm += -LAND_intercp;

  // // Evapo Interception
  NumericVector LAND_intercepEvapo_mm = evatransActual_SupplyPow(ATMOS_potentialEvatrans_mm,
                                                                 LAND_interceptWater_mm, LAND_interceptCapacity_mm,
                                                                 param_EVATRANS_sup_k, param_EVATRANS_sup_gamma);
  LAND_interceptWater_mm += -LAND_intercepEvapo_mm;
  ATMOS_potentialEvatrans_mm += -LAND_intercepEvapo_mm;

  // // Built up
  NumericVector LAND_runoffBuiltup_mm = ATMOS_rainFall_mm * (LAND_builtRatio_1 * .5);
  ATMOS_rainFall_mm += -LAND_runoffBuiltup_mm;


  // // Evapo soil
  SOIL_evatrans_mm = evatransActual_WaterGAP(ATMOS_potentialEvatrans_mm, SOIL_water_mm, SOIL_capacity_mm, param_EVATRANS_wat_petmax);
  SOIL_water_mm += -SOIL_evatrans_mm;
  LAND_water_mm = ATMOS_rainFall_mm;


  SOIL_evatrans_mm += LAND_interceptWater_mm;
  // // Snow melt
  SNOW_melt_mm = snowMelt_Factor(SNOW_ice_mm, ATMOS_temperature_Cel, param_SNOW_fac_f, param_SNOW_fac_Tmelt);
  LAND_water_mm += SNOW_melt_mm;
  SNOW_ice_mm += -SNOW_melt_mm;
  SNOW_ice_mm += ATMOS_snowFall_mm;

  // // soil infiltration
  SOIL_INFILT_mm = infilt_HBV(LAND_water_mm, SOIL_water_mm, SOIL_capacity_mm, param_INFILT_hbv_beta);
  SOIL_water_mm += SOIL_INFILT_mm;
  LAND_runoff_mm = LAND_water_mm - SOIL_INFILT_mm + LAND_runoffBuiltup_mm;

  // // soil percolation
  SOIL_percola_mm = percola_Arno(SOIL_water_mm, SOIL_capacity_mm, SOIL_potentialPercola_mm, param_PERCOLA_arn_thresh, param_PERCOLA_arn_k);
  GROUND_water_mm += SOIL_percola_mm;
  SOIL_water_mm += - SOIL_percola_mm;

  // // ground water
  GROUND_basefloW_mm = baseflow_SupplyRatio(GROUND_water_mm, param_BASEFLOW_sur_k);
  GROUND_water_mm += - GROUND_basefloW_mm;



  return (LAND_runoff_mm + GROUND_basefloW_mm) * CELL_landArea_km2 * 1000;
}
