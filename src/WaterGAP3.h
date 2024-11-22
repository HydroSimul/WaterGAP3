// Defines a header file containing function for WaterGAP3_HL/
#ifndef WATERGAP3_H
#define WATERGAP3_H

#include "00utilis.h"

NumericVector atmosSnow_ThresholdT(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector param_ATMOS_thr_Ts
);
NumericVector evatransPotential_TurcWendling(
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector param_EVATRANS_tur_k
);
NumericVector evatransActual_UBC(
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_EVATRANS_ubc_gamma
);
NumericVector evatransActual_VIC(
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_EVATRANS_vic_gamma
);
NumericVector snowMelt_Factor(
    NumericVector SNOW_ice_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector param_SNOW_fac_f,
    NumericVector param_SNOW_fac_Tmelt
);

NumericVector infilt_HBV(
    NumericVector LAND_water_mm,
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector param_INFILT_hbv_beta
);
NumericVector percola_Arno(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialPercola_mm,
    NumericVector param_PERCOLA_arn_thresh,
    NumericVector param_PERCOLA_arn_k
);
NumericVector baseflow_GR4Jfix(
    NumericVector GROUND_water_mm,
    NumericVector GROUND_capacity_mm,
    NumericVector param_BASEFLOW_grf_gamma
);

NumericVector lake_AcceptPow(
  NumericVector Lake_water_m3,
  NumericVector Lake_inflow_m3,
  NumericVector Lake_capacity_m3,
  NumericVector param_Lake_acp_storeFactor,
  NumericVector param_Lake_acp_gamma
);

NumericVector confluen_WaterGAP3(
    NumericVector CONFLUEN_cellInflow_m3,
    NumericVector &RIVER_water_m3,
    NumericVector RIVER_length_km,
    NumericVector RIVER_velocity_km,
    List celL_cellNumberStep_int,
    List celL_inflowCellNumberStep_int
);


NumericVector confluen_WaterGAP3_L(
    NumericVector CONFLUEN_cellInflow_m3,
    NumericVector &RIVER_water_m3,
    NumericVector RIVER_length_km,
    NumericVector RIVER_velocity_km,
    List celL_cellNumberStep_int,
    List celL_inflowCellNumberStep_int,
    IntegerVector Riverlak_cellNumber_int,
    NumericVector &Riverlak_water_m3,
    NumericVector Riverlak_capacity_m3,
    NumericVector param_Riverlak_lin_storeFactor
);

NumericVector module_land_WaterGAP3(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector& ATMOS_snowFall_mm,
    NumericVector& SNOW_ice_mm,
    NumericVector& LAND_runoff_mm,
    NumericVector& SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialPercola_mm,
    NumericVector& SOIL_EVATRANS_mm,
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
    NumericVector param_BASEFLOW_grf_gamma);

#endif

