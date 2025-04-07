#ifndef __MODEL__
#define __MODEL__

#include "00utilis.h"

NumericMatrix WaterGAP3_N(
    std::string name_Region,
    int n_time,
    int n_spat,
    NumericMatrix ATMOS_precipitation_mm,
    NumericMatrix ATMOS_temperature_Cel,
    NumericMatrix ATMOS_solarRadiat_MJ,
    NumericMatrix ATMOS_solarRadiatClearSky_MJ,
    IntegerVector Upstream_cellNumber_int,
    NumericMatrix Upstream_streamflow_m3,
    NumericVector SNOW_ice_mm,
    NumericVector LAND_area_km2,
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
    NumericVector param_Lake_acp_storeFactor,
    NumericVector param_Lake_acp_gamma,
    NumericVector param_Riverlak_lin_storeFactor,
    std::string path_FinalState = "NonExport",
    std::string path_VariExport = "NonExport"
);
#endif // __MODEL__
