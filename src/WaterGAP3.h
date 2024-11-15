// Defines a header file containing function for WaterGAP3_HL/
#ifndef WATERGAP3_H
#define WATERGAP3_H

#include "00utilis.h"

NumericVector atmosSnow_ThresholdT(
    NumericVector atmos_precipitation_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector param_atmos_thr_Ts
);
NumericVector evatransPotential_TurcWendling(
    NumericVector atmos_temperature_Cel,
    NumericVector atmos_solarRadiat_MJ,
    NumericVector param_evatrans_tur_k
);
NumericVector evatransActual_UBC(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_ubc_gamma
);
NumericVector evatransActual_VIC(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_vic_gamma
);
NumericVector snowMelt_Factor(
    NumericVector snow_ice_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector param_snow_fac_f,
    NumericVector param_snow_fac_Tmelt
);

NumericVector infilt_HBV(
    NumericVector land_water_mm,
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector param_infilt_hbv_beta
);
NumericVector percola_Arno(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_potentialPercola_mm,
    NumericVector param_percola_arn_thresh,
    NumericVector param_percola_arn_k
);
NumericVector baseflow_GR4Jfix(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_baseflow_grf_gamma
);
NumericVector confluen_IUH2S(
    NumericVector land_runoff_mm,
    NumericVector ground_baseflow_mm,
    NumericVector confluen_iuhLand_1,
    NumericVector confluen_iuhGround_1
);

NumericVector lake_AcceptPow(
  NumericVector lake_water_m3,
  NumericVector lake_inflow_m3,
  NumericVector lake_capacity_m3,
  NumericVector param_lake_acp_storeFactor,
  NumericVector param_lake_acp_gamma
);


NumericVector confluen_WaterGAP3_L(
    NumericVector confluen_cellInflow_m3,
    NumericVector &river_water_m3,
    NumericVector river_length_km,
    NumericVector river_velocity_km,
    IntegerVector riverlake_cellNumber_int,
    NumericVector &riverlake_water_m3,
    NumericVector riverlake_capacity_m3,
    List basin_cellNumberStep_int,
    List basin_inflowCellNumberStep_int,
    NumericVector param_riverlake_lin_storeFactor
);


NumericVector confluen_WaterGAP3_LR(
    NumericVector confluen_cellInflow_m3,
    NumericVector &river_water_m3,
    NumericVector river_length_km,
    NumericVector river_velocity_km,
    IntegerVector riverlake_cellNumber_int,
    NumericVector &riverlake_water_m3,
    NumericVector riverlake_capacity_m3,
    IntegerVector reservoir_cellNumber_int,
    NumericVector &reservoir_water_m3,
    NumericVector reservoir_capacity_m3,
    NumericVector reservoir_demand_m3,
    NumericVector reservoir_yearInflow_m3,
    NumericVector reservoir_yearDemand_m3,
    NumericVector &reservoir_yearRelase_m3,
    LogicalVector reservoir_isOperateStart_01,
    LogicalVector reservoir_isIrrigate_01,
    List basin_cellNumberStep_int,
    List basin_inflowCellNumberStep_int,
    NumericVector param_riverlake_lin_storeFactor,
    NumericVector param_reservoir_han_alpha,
    NumericVector param_reservoir_han_kDemand
);
#endif
