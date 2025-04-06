// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/WaterGAP3.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// WaterGAP3_N
List WaterGAP3_N(int n_time, int n_spat, NumericMatrix ATMOS_precipitation_mm, NumericMatrix ATMOS_temperature_Cel, NumericMatrix ATMOS_solarRadiat_MJ, NumericMatrix ATMOS_solarRadiatClearSky_MJ, IntegerVector Upstream_cellNumber_int, NumericMatrix Upstream_streamflow_m3, NumericVector SNOW_ice_mm, NumericVector LAND_area_km2, NumericVector LAND_albedo_1, NumericVector LAND_snowAlbedo_1, NumericVector LAND_builtRatio_1, NumericVector LAND_interceptWater_mm, NumericMatrix LAND_interceptCapacity_mm, NumericVector SOIL_water_mm, NumericVector SOIL_capacity_mm, NumericVector SOIL_potentialPercola_mm, NumericVector GROUND_water_mm, NumericVector RIVER_water_m3, NumericVector RIVER_length_km, NumericVector RIVER_velocity_km, NumericVector CELL_elevation_m, List CELL_cellNumberStep_int, List CELL_inflowCellNumberStep_int, IntegerVector Lake_cellNumber_int, NumericVector Lake_water_m3, NumericVector Lake_area_km2, NumericVector Lake_capacity_m3, NumericVector Lake_albedo_1, IntegerVector Riverlak_cellNumber_int, NumericVector Riverlak_water_m3, NumericVector Riverlak_area_km2, NumericVector Riverlak_capacity_m3, NumericVector Riverlak_albedo_1, NumericVector param_ATMOS_thr_Ts, NumericVector param_SNOW_fac_f, NumericVector param_SNOW_fac_Tmelt, NumericVector param_EVATRANS_prt_alpha, NumericVector param_EVATRANS_vic_gamma, NumericVector param_EVATRANS_sup_k, NumericVector param_EVATRANS_sup_gamma, NumericVector param_EVATRANS_wat_petmax, NumericVector param_INFILT_hbv_beta, LogicalVector param_PERCOLA_wat_01, NumericVector param_PERCOLA_wat_k, NumericVector param_PERCOLA_wat_thresh, NumericVector param_BASEFLOW_sur_k, NumericVector param_Lake_Eva_vic_gamma, NumericVector param_Lake_acp_storeFactor, NumericVector param_Lake_acp_gamma, NumericVector param_Riverlak_Eva_vic_gamma, NumericVector param_Riverlak_lin_storeFactor, bool if_allVariExport);
static SEXP _WaterGAP3_WaterGAP3_N_try(SEXP n_timeSEXP, SEXP n_spatSEXP, SEXP ATMOS_precipitation_mmSEXP, SEXP ATMOS_temperature_CelSEXP, SEXP ATMOS_solarRadiat_MJSEXP, SEXP ATMOS_solarRadiatClearSky_MJSEXP, SEXP Upstream_cellNumber_intSEXP, SEXP Upstream_streamflow_m3SEXP, SEXP SNOW_ice_mmSEXP, SEXP LAND_area_km2SEXP, SEXP LAND_albedo_1SEXP, SEXP LAND_snowAlbedo_1SEXP, SEXP LAND_builtRatio_1SEXP, SEXP LAND_interceptWater_mmSEXP, SEXP LAND_interceptCapacity_mmSEXP, SEXP SOIL_water_mmSEXP, SEXP SOIL_capacity_mmSEXP, SEXP SOIL_potentialPercola_mmSEXP, SEXP GROUND_water_mmSEXP, SEXP RIVER_water_m3SEXP, SEXP RIVER_length_kmSEXP, SEXP RIVER_velocity_kmSEXP, SEXP CELL_elevation_mSEXP, SEXP CELL_cellNumberStep_intSEXP, SEXP CELL_inflowCellNumberStep_intSEXP, SEXP Lake_cellNumber_intSEXP, SEXP Lake_water_m3SEXP, SEXP Lake_area_km2SEXP, SEXP Lake_capacity_m3SEXP, SEXP Lake_albedo_1SEXP, SEXP Riverlak_cellNumber_intSEXP, SEXP Riverlak_water_m3SEXP, SEXP Riverlak_area_km2SEXP, SEXP Riverlak_capacity_m3SEXP, SEXP Riverlak_albedo_1SEXP, SEXP param_ATMOS_thr_TsSEXP, SEXP param_SNOW_fac_fSEXP, SEXP param_SNOW_fac_TmeltSEXP, SEXP param_EVATRANS_prt_alphaSEXP, SEXP param_EVATRANS_vic_gammaSEXP, SEXP param_EVATRANS_sup_kSEXP, SEXP param_EVATRANS_sup_gammaSEXP, SEXP param_EVATRANS_wat_petmaxSEXP, SEXP param_INFILT_hbv_betaSEXP, SEXP param_PERCOLA_wat_01SEXP, SEXP param_PERCOLA_wat_kSEXP, SEXP param_PERCOLA_wat_threshSEXP, SEXP param_BASEFLOW_sur_kSEXP, SEXP param_Lake_Eva_vic_gammaSEXP, SEXP param_Lake_acp_storeFactorSEXP, SEXP param_Lake_acp_gammaSEXP, SEXP param_Riverlak_Eva_vic_gammaSEXP, SEXP param_Riverlak_lin_storeFactorSEXP, SEXP if_allVariExportSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n_time(n_timeSEXP);
    Rcpp::traits::input_parameter< int >::type n_spat(n_spatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ATMOS_precipitation_mm(ATMOS_precipitation_mmSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ATMOS_temperature_Cel(ATMOS_temperature_CelSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ATMOS_solarRadiat_MJ(ATMOS_solarRadiat_MJSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ATMOS_solarRadiatClearSky_MJ(ATMOS_solarRadiatClearSky_MJSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Upstream_cellNumber_int(Upstream_cellNumber_intSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Upstream_streamflow_m3(Upstream_streamflow_m3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SNOW_ice_mm(SNOW_ice_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LAND_area_km2(LAND_area_km2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LAND_albedo_1(LAND_albedo_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LAND_snowAlbedo_1(LAND_snowAlbedo_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LAND_builtRatio_1(LAND_builtRatio_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LAND_interceptWater_mm(LAND_interceptWater_mmSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LAND_interceptCapacity_mm(LAND_interceptCapacity_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SOIL_water_mm(SOIL_water_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SOIL_capacity_mm(SOIL_capacity_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SOIL_potentialPercola_mm(SOIL_potentialPercola_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type GROUND_water_mm(GROUND_water_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type RIVER_water_m3(RIVER_water_m3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type RIVER_length_km(RIVER_length_kmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type RIVER_velocity_km(RIVER_velocity_kmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type CELL_elevation_m(CELL_elevation_mSEXP);
    Rcpp::traits::input_parameter< List >::type CELL_cellNumberStep_int(CELL_cellNumberStep_intSEXP);
    Rcpp::traits::input_parameter< List >::type CELL_inflowCellNumberStep_int(CELL_inflowCellNumberStep_intSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Lake_cellNumber_int(Lake_cellNumber_intSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lake_water_m3(Lake_water_m3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lake_area_km2(Lake_area_km2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lake_capacity_m3(Lake_capacity_m3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lake_albedo_1(Lake_albedo_1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Riverlak_cellNumber_int(Riverlak_cellNumber_intSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Riverlak_water_m3(Riverlak_water_m3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Riverlak_area_km2(Riverlak_area_km2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Riverlak_capacity_m3(Riverlak_capacity_m3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Riverlak_albedo_1(Riverlak_albedo_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_ATMOS_thr_Ts(param_ATMOS_thr_TsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_SNOW_fac_f(param_SNOW_fac_fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_SNOW_fac_Tmelt(param_SNOW_fac_TmeltSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_EVATRANS_prt_alpha(param_EVATRANS_prt_alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_EVATRANS_vic_gamma(param_EVATRANS_vic_gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_EVATRANS_sup_k(param_EVATRANS_sup_kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_EVATRANS_sup_gamma(param_EVATRANS_sup_gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_EVATRANS_wat_petmax(param_EVATRANS_wat_petmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_INFILT_hbv_beta(param_INFILT_hbv_betaSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type param_PERCOLA_wat_01(param_PERCOLA_wat_01SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_PERCOLA_wat_k(param_PERCOLA_wat_kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_PERCOLA_wat_thresh(param_PERCOLA_wat_threshSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_BASEFLOW_sur_k(param_BASEFLOW_sur_kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_Lake_Eva_vic_gamma(param_Lake_Eva_vic_gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_Lake_acp_storeFactor(param_Lake_acp_storeFactorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_Lake_acp_gamma(param_Lake_acp_gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_Riverlak_Eva_vic_gamma(param_Riverlak_Eva_vic_gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_Riverlak_lin_storeFactor(param_Riverlak_lin_storeFactorSEXP);
    Rcpp::traits::input_parameter< bool >::type if_allVariExport(if_allVariExportSEXP);
    rcpp_result_gen = Rcpp::wrap(WaterGAP3_N(n_time, n_spat, ATMOS_precipitation_mm, ATMOS_temperature_Cel, ATMOS_solarRadiat_MJ, ATMOS_solarRadiatClearSky_MJ, Upstream_cellNumber_int, Upstream_streamflow_m3, SNOW_ice_mm, LAND_area_km2, LAND_albedo_1, LAND_snowAlbedo_1, LAND_builtRatio_1, LAND_interceptWater_mm, LAND_interceptCapacity_mm, SOIL_water_mm, SOIL_capacity_mm, SOIL_potentialPercola_mm, GROUND_water_mm, RIVER_water_m3, RIVER_length_km, RIVER_velocity_km, CELL_elevation_m, CELL_cellNumberStep_int, CELL_inflowCellNumberStep_int, Lake_cellNumber_int, Lake_water_m3, Lake_area_km2, Lake_capacity_m3, Lake_albedo_1, Riverlak_cellNumber_int, Riverlak_water_m3, Riverlak_area_km2, Riverlak_capacity_m3, Riverlak_albedo_1, param_ATMOS_thr_Ts, param_SNOW_fac_f, param_SNOW_fac_Tmelt, param_EVATRANS_prt_alpha, param_EVATRANS_vic_gamma, param_EVATRANS_sup_k, param_EVATRANS_sup_gamma, param_EVATRANS_wat_petmax, param_INFILT_hbv_beta, param_PERCOLA_wat_01, param_PERCOLA_wat_k, param_PERCOLA_wat_thresh, param_BASEFLOW_sur_k, param_Lake_Eva_vic_gamma, param_Lake_acp_storeFactor, param_Lake_acp_gamma, param_Riverlak_Eva_vic_gamma, param_Riverlak_lin_storeFactor, if_allVariExport));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_WaterGAP3_N(SEXP n_timeSEXP, SEXP n_spatSEXP, SEXP ATMOS_precipitation_mmSEXP, SEXP ATMOS_temperature_CelSEXP, SEXP ATMOS_solarRadiat_MJSEXP, SEXP ATMOS_solarRadiatClearSky_MJSEXP, SEXP Upstream_cellNumber_intSEXP, SEXP Upstream_streamflow_m3SEXP, SEXP SNOW_ice_mmSEXP, SEXP LAND_area_km2SEXP, SEXP LAND_albedo_1SEXP, SEXP LAND_snowAlbedo_1SEXP, SEXP LAND_builtRatio_1SEXP, SEXP LAND_interceptWater_mmSEXP, SEXP LAND_interceptCapacity_mmSEXP, SEXP SOIL_water_mmSEXP, SEXP SOIL_capacity_mmSEXP, SEXP SOIL_potentialPercola_mmSEXP, SEXP GROUND_water_mmSEXP, SEXP RIVER_water_m3SEXP, SEXP RIVER_length_kmSEXP, SEXP RIVER_velocity_kmSEXP, SEXP CELL_elevation_mSEXP, SEXP CELL_cellNumberStep_intSEXP, SEXP CELL_inflowCellNumberStep_intSEXP, SEXP Lake_cellNumber_intSEXP, SEXP Lake_water_m3SEXP, SEXP Lake_area_km2SEXP, SEXP Lake_capacity_m3SEXP, SEXP Lake_albedo_1SEXP, SEXP Riverlak_cellNumber_intSEXP, SEXP Riverlak_water_m3SEXP, SEXP Riverlak_area_km2SEXP, SEXP Riverlak_capacity_m3SEXP, SEXP Riverlak_albedo_1SEXP, SEXP param_ATMOS_thr_TsSEXP, SEXP param_SNOW_fac_fSEXP, SEXP param_SNOW_fac_TmeltSEXP, SEXP param_EVATRANS_prt_alphaSEXP, SEXP param_EVATRANS_vic_gammaSEXP, SEXP param_EVATRANS_sup_kSEXP, SEXP param_EVATRANS_sup_gammaSEXP, SEXP param_EVATRANS_wat_petmaxSEXP, SEXP param_INFILT_hbv_betaSEXP, SEXP param_PERCOLA_wat_01SEXP, SEXP param_PERCOLA_wat_kSEXP, SEXP param_PERCOLA_wat_threshSEXP, SEXP param_BASEFLOW_sur_kSEXP, SEXP param_Lake_Eva_vic_gammaSEXP, SEXP param_Lake_acp_storeFactorSEXP, SEXP param_Lake_acp_gammaSEXP, SEXP param_Riverlak_Eva_vic_gammaSEXP, SEXP param_Riverlak_lin_storeFactorSEXP, SEXP if_allVariExportSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_WaterGAP3_N_try(n_timeSEXP, n_spatSEXP, ATMOS_precipitation_mmSEXP, ATMOS_temperature_CelSEXP, ATMOS_solarRadiat_MJSEXP, ATMOS_solarRadiatClearSky_MJSEXP, Upstream_cellNumber_intSEXP, Upstream_streamflow_m3SEXP, SNOW_ice_mmSEXP, LAND_area_km2SEXP, LAND_albedo_1SEXP, LAND_snowAlbedo_1SEXP, LAND_builtRatio_1SEXP, LAND_interceptWater_mmSEXP, LAND_interceptCapacity_mmSEXP, SOIL_water_mmSEXP, SOIL_capacity_mmSEXP, SOIL_potentialPercola_mmSEXP, GROUND_water_mmSEXP, RIVER_water_m3SEXP, RIVER_length_kmSEXP, RIVER_velocity_kmSEXP, CELL_elevation_mSEXP, CELL_cellNumberStep_intSEXP, CELL_inflowCellNumberStep_intSEXP, Lake_cellNumber_intSEXP, Lake_water_m3SEXP, Lake_area_km2SEXP, Lake_capacity_m3SEXP, Lake_albedo_1SEXP, Riverlak_cellNumber_intSEXP, Riverlak_water_m3SEXP, Riverlak_area_km2SEXP, Riverlak_capacity_m3SEXP, Riverlak_albedo_1SEXP, param_ATMOS_thr_TsSEXP, param_SNOW_fac_fSEXP, param_SNOW_fac_TmeltSEXP, param_EVATRANS_prt_alphaSEXP, param_EVATRANS_vic_gammaSEXP, param_EVATRANS_sup_kSEXP, param_EVATRANS_sup_gammaSEXP, param_EVATRANS_wat_petmaxSEXP, param_INFILT_hbv_betaSEXP, param_PERCOLA_wat_01SEXP, param_PERCOLA_wat_kSEXP, param_PERCOLA_wat_threshSEXP, param_BASEFLOW_sur_kSEXP, param_Lake_Eva_vic_gammaSEXP, param_Lake_acp_storeFactorSEXP, param_Lake_acp_gammaSEXP, param_Riverlak_Eva_vic_gammaSEXP, param_Riverlak_lin_storeFactorSEXP, if_allVariExportSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// write_nc_WG3
void write_nc_WG3(NumericMatrix mat_Data_WG3, IntegerVector dim_Time, IntegerVector dim_Spat, std::string path_File, std::string name_Variable, std::string str_Continent, std::string suffix_File);
static SEXP _WaterGAP3_write_nc_WG3_try(SEXP mat_Data_WG3SEXP, SEXP dim_TimeSEXP, SEXP dim_SpatSEXP, SEXP path_FileSEXP, SEXP name_VariableSEXP, SEXP str_ContinentSEXP, SEXP suffix_FileSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_Data_WG3(mat_Data_WG3SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dim_Time(dim_TimeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dim_Spat(dim_SpatSEXP);
    Rcpp::traits::input_parameter< std::string >::type path_File(path_FileSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_Variable(name_VariableSEXP);
    Rcpp::traits::input_parameter< std::string >::type str_Continent(str_ContinentSEXP);
    Rcpp::traits::input_parameter< std::string >::type suffix_File(suffix_FileSEXP);
    write_nc_WG3(mat_Data_WG3, dim_Time, dim_Spat, path_File, name_Variable, str_Continent, suffix_File);
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_write_nc_WG3(SEXP mat_Data_WG3SEXP, SEXP dim_TimeSEXP, SEXP dim_SpatSEXP, SEXP path_FileSEXP, SEXP name_VariableSEXP, SEXP str_ContinentSEXP, SEXP suffix_FileSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_write_nc_WG3_try(mat_Data_WG3SEXP, dim_TimeSEXP, dim_SpatSEXP, path_FileSEXP, name_VariableSEXP, str_ContinentSEXP, suffix_FileSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// read_nc_WG3
NumericMatrix read_nc_WG3(std::string path_File, std::string name_Variable);
static SEXP _WaterGAP3_read_nc_WG3_try(SEXP path_FileSEXP, SEXP name_VariableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type path_File(path_FileSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_Variable(name_VariableSEXP);
    rcpp_result_gen = Rcpp::wrap(read_nc_WG3(path_File, name_Variable));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_read_nc_WG3(SEXP path_FileSEXP, SEXP name_VariableSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_read_nc_WG3_try(path_FileSEXP, name_VariableSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bind_nc_WG3
void bind_nc_WG3(std::vector<std::string> path_FilesBind, std::string path_Out, std::string name_Variable);
static SEXP _WaterGAP3_bind_nc_WG3_try(SEXP path_FilesBindSEXP, SEXP path_OutSEXP, SEXP name_VariableSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< std::vector<std::string> >::type path_FilesBind(path_FilesBindSEXP);
    Rcpp::traits::input_parameter< std::string >::type path_Out(path_OutSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_Variable(name_VariableSEXP);
    bind_nc_WG3(path_FilesBind, path_Out, name_Variable);
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_bind_nc_WG3(SEXP path_FilesBindSEXP, SEXP path_OutSEXP, SEXP name_VariableSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_bind_nc_WG3_try(path_FilesBindSEXP, path_OutSEXP, name_VariableSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// read_nc_dim_WG3
Rcpp::IntegerVector read_nc_dim_WG3(std::string path_File, std::string dim_name);
static SEXP _WaterGAP3_read_nc_dim_WG3_try(SEXP path_FileSEXP, SEXP dim_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type path_File(path_FileSEXP);
    Rcpp::traits::input_parameter< std::string >::type dim_name(dim_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_nc_dim_WG3(path_File, dim_name));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_read_nc_dim_WG3(SEXP path_FileSEXP, SEXP dim_nameSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_read_nc_dim_WG3_try(path_FileSEXP, dim_nameSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// read_unf
SEXP read_unf(std::string fn_UNF);
static SEXP _WaterGAP3_read_unf_try(SEXP fn_UNFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type fn_UNF(fn_UNFSEXP);
    rcpp_result_gen = Rcpp::wrap(read_unf(fn_UNF));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_read_unf(SEXP fn_UNFSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_read_unf_try(fn_UNFSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// write_unf
void write_unf(SEXP data_Export, std::string fn_UNF);
static SEXP _WaterGAP3_write_unf_try(SEXP data_ExportSEXP, SEXP fn_UNFSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< SEXP >::type data_Export(data_ExportSEXP);
    Rcpp::traits::input_parameter< std::string >::type fn_UNF(fn_UNFSEXP);
    write_unf(data_Export, fn_UNF);
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_write_unf(SEXP data_ExportSEXP, SEXP fn_UNFSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_write_unf_try(data_ExportSEXP, fn_UNFSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// save_wgmat
void save_wgmat(const NumericMatrix& matrix, const std::string& filename);
static SEXP _WaterGAP3_save_wgmat_try(SEXP matrixSEXP, SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< const NumericMatrix& >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    save_wgmat(matrix, filename);
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_save_wgmat(SEXP matrixSEXP, SEXP filenameSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_save_wgmat_try(matrixSEXP, filenameSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// load_wgmat
NumericMatrix load_wgmat(const std::string& filename);
static SEXP _WaterGAP3_load_wgmat_try(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(load_wgmat(filename));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_load_wgmat(SEXP filenameSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_load_wgmat_try(filenameSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bind_wgmat
void bind_wgmat(const StringVector& input_files, const std::string& output_file);
static SEXP _WaterGAP3_bind_wgmat_try(SEXP input_filesSEXP, SEXP output_fileSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< const StringVector& >::type input_files(input_filesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type output_file(output_fileSEXP);
    bind_wgmat(input_files, output_file);
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _WaterGAP3_bind_wgmat(SEXP input_filesSEXP, SEXP output_fileSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_WaterGAP3_bind_wgmat_try(input_filesSEXP, output_fileSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _WaterGAP3_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("List(*WaterGAP3_N)(int,int,NumericMatrix,NumericMatrix,NumericMatrix,NumericMatrix,IntegerVector,NumericMatrix,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericMatrix,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,List,List,IntegerVector,NumericVector,NumericVector,NumericVector,NumericVector,IntegerVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,LogicalVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,bool)");
        signatures.insert("void(*write_nc_WG3)(NumericMatrix,IntegerVector,IntegerVector,std::string,std::string,std::string,std::string)");
        signatures.insert("NumericMatrix(*read_nc_WG3)(std::string,std::string)");
        signatures.insert("void(*bind_nc_WG3)(std::vector<std::string>,std::string,std::string)");
        signatures.insert("Rcpp::IntegerVector(*read_nc_dim_WG3)(std::string,std::string)");
        signatures.insert("SEXP(*read_unf)(std::string)");
        signatures.insert("void(*write_unf)(SEXP,std::string)");
        signatures.insert("void(*save_wgmat)(const NumericMatrix&,const std::string&)");
        signatures.insert("NumericMatrix(*load_wgmat)(const std::string&)");
        signatures.insert("void(*bind_wgmat)(const StringVector&,const std::string&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _WaterGAP3_RcppExport_registerCCallable() { 
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_WaterGAP3_N", (DL_FUNC)_WaterGAP3_WaterGAP3_N_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_write_nc_WG3", (DL_FUNC)_WaterGAP3_write_nc_WG3_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_read_nc_WG3", (DL_FUNC)_WaterGAP3_read_nc_WG3_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_bind_nc_WG3", (DL_FUNC)_WaterGAP3_bind_nc_WG3_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_read_nc_dim_WG3", (DL_FUNC)_WaterGAP3_read_nc_dim_WG3_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_read_unf", (DL_FUNC)_WaterGAP3_read_unf_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_write_unf", (DL_FUNC)_WaterGAP3_write_unf_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_save_wgmat", (DL_FUNC)_WaterGAP3_save_wgmat_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_load_wgmat", (DL_FUNC)_WaterGAP3_load_wgmat_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_bind_wgmat", (DL_FUNC)_WaterGAP3_bind_wgmat_try);
    R_RegisterCCallable("WaterGAP3", "_WaterGAP3_RcppExport_validate", (DL_FUNC)_WaterGAP3_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_WaterGAP3_WaterGAP3_N", (DL_FUNC) &_WaterGAP3_WaterGAP3_N, 54},
    {"_WaterGAP3_write_nc_WG3", (DL_FUNC) &_WaterGAP3_write_nc_WG3, 7},
    {"_WaterGAP3_read_nc_WG3", (DL_FUNC) &_WaterGAP3_read_nc_WG3, 2},
    {"_WaterGAP3_bind_nc_WG3", (DL_FUNC) &_WaterGAP3_bind_nc_WG3, 3},
    {"_WaterGAP3_read_nc_dim_WG3", (DL_FUNC) &_WaterGAP3_read_nc_dim_WG3, 2},
    {"_WaterGAP3_read_unf", (DL_FUNC) &_WaterGAP3_read_unf, 1},
    {"_WaterGAP3_write_unf", (DL_FUNC) &_WaterGAP3_write_unf, 2},
    {"_WaterGAP3_save_wgmat", (DL_FUNC) &_WaterGAP3_save_wgmat, 2},
    {"_WaterGAP3_load_wgmat", (DL_FUNC) &_WaterGAP3_load_wgmat, 1},
    {"_WaterGAP3_bind_wgmat", (DL_FUNC) &_WaterGAP3_bind_wgmat, 2},
    {"_WaterGAP3_RcppExport_registerCCallable", (DL_FUNC) &_WaterGAP3_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_WaterGAP3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
