// Defines a header file containing function for WATERGAP3_Process/
#ifndef WATERGAP3_Module
#define WATERGAP3_Module

#include "00utilis.h"

NumericVector module_land_WaterGAP3(
    NumericVector ATMOS_precipitation_mm,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector ATMOS_solarRadiatClearSky_MJ,
    NumericVector& ATMOS_potentialEvatrans_mm,
    NumericVector& ATMOS_snowFall_mm,
    NumericVector& SNOW_ice_mm,
    NumericVector LAND_area_km2,
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
    NumericVector param_BASEFLOW_sur_k);

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
    NumericVector param_EVATRANS_vic_gamma);

NumericVector module_lake_WaterGAP3(
    NumericVector& Lake_water_m3,
    NumericVector Lake_capacity_m3,
    NumericVector Lake_verticalInflow_m3,
    NumericVector Lake_inflow_m3,
    NumericVector param_Lake_acp_storeFactor,
    NumericVector param_Lake_acp_gamma);
#endif

