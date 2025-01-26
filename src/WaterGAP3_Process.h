// Defines a header file containing function for WATERGAP3_Process/
#ifndef WATERGAP3_Process
#define WATERGAP3_Process

#include "00utilis.h"

NumericVector atmosSnow_ThresholdT(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector param_ATMOS_thr_Ts
);
NumericVector intercep_Full(
    NumericVector ATMOS_precipitation_mm,
    NumericVector LAND_interceptWater_mm,
    NumericVector LAND_interceptCapacity_mm
);
NumericVector evatransActual_SupplyPow(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_EVATRANS_sup_k,
    NumericVector param_EVATRANS_sup_gamma
);
NumericVector evatransActual_SupplyRatio(
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_EVATRANS_sur_k
);
NumericVector evatransActual_WaterGAP3(
    NumericVector ATMOS_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_EVATRANS_wat_petmax
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
NumericVector percola_WaterGAP3(
    NumericVector Land_water_mm,
    NumericVector SOIL_potentialPercola_mm,
    LogicalVector param_PERCOLA_wat_01,
    NumericVector param_PERCOLA_wat_thresh,
    NumericVector param_PERCOLA_wat_k
);
NumericVector baseflow_GR4Jfix(
    NumericVector GROUND_water_mm,
    NumericVector GROUND_capacity_mm,
    NumericVector param_BASEFLOW_grf_gamma
);
NumericVector baseflow_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector param_baseflow_sur_k
);
NumericVector lake_AcceptPow(
    NumericVector Lake_water_m3,
    NumericVector Lake_capacity_m3,
    NumericVector param_Lake_acp_storeFactor,
    NumericVector param_Lake_acp_gamma
);

NumericVector confluen_WaterGAP3(
    NumericVector &RIVER_water_m3,
    NumericVector RIVER_length_km,
    NumericVector RIVER_velocity_km,
    NumericVector RIVER_outflow_m3,
    List CELL_cellNumberStep_int,
    List CELL_inflowCellNumberStep_int
);


NumericVector confluen_WaterGAP3_L(
    NumericVector &RIVER_water_m3,
    NumericVector RIVER_length_km,
    NumericVector RIVER_velocity_km,
    NumericVector RIVER_outflow_m3,
    List CELL_cellNumberStep_int,
    List CELL_inflowCellNumberStep_int,
    IntegerVector Riverlak_cellNumber_int,
    NumericVector &Riverlak_water_m3,
    NumericVector Riverlak_capacity_m3,
    NumericVector param_Riverlak_lin_storeFactor
);

NumericVector confluen_WaterGAP3_LR(
    NumericVector CONFLUEN_cellInflow_m3,
    NumericVector &RIVER_water_m3,
    NumericVector RIVER_length_km,
    NumericVector RIVER_velocity_km,
    IntegerVector Riverlak_cellNumber_int,
    NumericVector &Riverlak_water_m3,
    NumericVector Riverlak_capacity_m3,
    IntegerVector Reservoi_cellNumber_int,
    NumericVector &Reservoi_water_m3,
    NumericVector Reservoi_capacity_m3,
    NumericVector Reservoi_demand_m3,
    NumericVector Reservoi_yearInflow_m3,
    NumericVector Reservoi_yearDemand_m3,
    NumericVector &Reservoi_yearRelase_m3,
    LogicalVector Reservoi_isOperateStart_01,
    LogicalVector Reservoi_isIrrigate_01,
    List CELL_cellNumberStep_int,
    List CELL_inflowCellNumberStep_int,
    NumericVector param_Riverlak_lin_storeFactor,
    NumericVector param_Reservoi_han_alpha,
    NumericVector param_Reservoi_han_kDemand
);

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
    LogicalVector param_PERCOLA_wat_01,
    NumericVector param_PERCOLA_wat_thresh,
    NumericVector param_PERCOLA_wat_k,
    NumericVector param_BASEFLOW_sur_k);

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

