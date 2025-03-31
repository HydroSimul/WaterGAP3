// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_WaterGAP3_RCPPEXPORTS_H_GEN_
#define RCPP_WaterGAP3_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace WaterGAP3 {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("WaterGAP3", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("WaterGAP3", "_WaterGAP3_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in WaterGAP3");
            }
        }
    }

    inline List WaterGAP3_H(int n_time, int n_spat, NumericMatrix ATMOS_precipitation_mm, NumericMatrix ATMOS_temperature_Cel, NumericMatrix ATMOS_potentialEvatrans_mm, IntegerVector Upstream_cellNumber_int, NumericMatrix Upstream_streamflow_m3, NumericVector SNOW_ice_mm, NumericVector LAND_builtRatio_1, NumericVector LAND_interceptWater_mm, NumericMatrix LAND_interceptCapacity_mm, NumericVector SOIL_water_mm, NumericVector SOIL_capacity_mm, NumericVector SOIL_potentialPercola_mm, NumericVector GROUND_water_mm, NumericVector RIVER_water_m3, NumericVector RIVER_length_km, NumericVector RIVER_velocity_km, NumericVector CELL_landArea_km2, List CELL_cellNumberStep_int, List CELL_inflowCellNumberStep_int, NumericVector param_ATMOS_thr_Ts, NumericVector param_SNOW_fac_f, NumericVector param_SNOW_fac_Tmelt, NumericVector param_EVATRANS_sup_k, NumericVector param_EVATRANS_sup_gamma, NumericVector param_EVATRANS_wat_petmax, NumericVector param_INFILT_hbv_beta, LogicalVector param_PERCOLA_wat_01, NumericVector param_PERCOLA_wat_k, NumericVector param_PERCOLA_wat_thresh, NumericVector param_BASEFLOW_sur_k, bool if_allVariExport = false) {
        typedef SEXP(*Ptr_WaterGAP3_H)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_WaterGAP3_H p_WaterGAP3_H = NULL;
        if (p_WaterGAP3_H == NULL) {
            validateSignature("List(*WaterGAP3_H)(int,int,NumericMatrix,NumericMatrix,NumericMatrix,IntegerVector,NumericMatrix,NumericVector,NumericVector,NumericVector,NumericMatrix,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,List,List,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,LogicalVector,NumericVector,NumericVector,NumericVector,bool)");
            p_WaterGAP3_H = (Ptr_WaterGAP3_H)R_GetCCallable("WaterGAP3", "_WaterGAP3_WaterGAP3_H");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_WaterGAP3_H(Shield<SEXP>(Rcpp::wrap(n_time)), Shield<SEXP>(Rcpp::wrap(n_spat)), Shield<SEXP>(Rcpp::wrap(ATMOS_precipitation_mm)), Shield<SEXP>(Rcpp::wrap(ATMOS_temperature_Cel)), Shield<SEXP>(Rcpp::wrap(ATMOS_potentialEvatrans_mm)), Shield<SEXP>(Rcpp::wrap(Upstream_cellNumber_int)), Shield<SEXP>(Rcpp::wrap(Upstream_streamflow_m3)), Shield<SEXP>(Rcpp::wrap(SNOW_ice_mm)), Shield<SEXP>(Rcpp::wrap(LAND_builtRatio_1)), Shield<SEXP>(Rcpp::wrap(LAND_interceptWater_mm)), Shield<SEXP>(Rcpp::wrap(LAND_interceptCapacity_mm)), Shield<SEXP>(Rcpp::wrap(SOIL_water_mm)), Shield<SEXP>(Rcpp::wrap(SOIL_capacity_mm)), Shield<SEXP>(Rcpp::wrap(SOIL_potentialPercola_mm)), Shield<SEXP>(Rcpp::wrap(GROUND_water_mm)), Shield<SEXP>(Rcpp::wrap(RIVER_water_m3)), Shield<SEXP>(Rcpp::wrap(RIVER_length_km)), Shield<SEXP>(Rcpp::wrap(RIVER_velocity_km)), Shield<SEXP>(Rcpp::wrap(CELL_landArea_km2)), Shield<SEXP>(Rcpp::wrap(CELL_cellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(CELL_inflowCellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(param_ATMOS_thr_Ts)), Shield<SEXP>(Rcpp::wrap(param_SNOW_fac_f)), Shield<SEXP>(Rcpp::wrap(param_SNOW_fac_Tmelt)), Shield<SEXP>(Rcpp::wrap(param_EVATRANS_sup_k)), Shield<SEXP>(Rcpp::wrap(param_EVATRANS_sup_gamma)), Shield<SEXP>(Rcpp::wrap(param_EVATRANS_wat_petmax)), Shield<SEXP>(Rcpp::wrap(param_INFILT_hbv_beta)), Shield<SEXP>(Rcpp::wrap(param_PERCOLA_wat_01)), Shield<SEXP>(Rcpp::wrap(param_PERCOLA_wat_k)), Shield<SEXP>(Rcpp::wrap(param_PERCOLA_wat_thresh)), Shield<SEXP>(Rcpp::wrap(param_BASEFLOW_sur_k)), Shield<SEXP>(Rcpp::wrap(if_allVariExport)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List WaterGAP3_HL(int n_time, int n_spat, NumericMatrix ATMOS_precipitation_mm, NumericMatrix ATMOS_temperature_Cel, NumericMatrix ATMOS_potentialEvatrans_mm, IntegerVector Upstream_cellNumber_int, NumericMatrix Upstream_streamflow_m3, NumericVector SNOW_ice_mm, NumericVector LAND_builtRatio_1, NumericVector LAND_interceptWater_mm, NumericMatrix LAND_interceptCapacity_mm, NumericVector SOIL_water_mm, NumericVector SOIL_capacity_mm, NumericVector SOIL_potentialPercola_mm, NumericVector GROUND_water_mm, NumericVector RIVER_water_m3, NumericVector RIVER_length_km, NumericVector RIVER_velocity_km, NumericVector CELL_landArea_km2, List CELL_cellNumberStep_int, List CELL_inflowCellNumberStep_int, IntegerVector Lake_cellNumber_int, NumericVector Lake_water_m3, NumericVector Lake_area_km2, NumericVector Lake_capacity_m3, IntegerVector Riverlak_cellNumber_int, NumericVector Riverlak_water_m3, NumericVector Riverlak_area_km2, NumericVector Riverlak_capacity_m3, NumericVector param_ATMOS_thr_Ts, NumericVector param_SNOW_fac_f, NumericVector param_SNOW_fac_Tmelt, NumericVector param_EVATRANS_sup_k, NumericVector param_EVATRANS_sup_gamma, NumericVector param_EVATRANS_wat_petmax, NumericVector param_INFILT_hbv_beta, LogicalVector param_PERCOLA_wat_01, NumericVector param_PERCOLA_wat_k, NumericVector param_PERCOLA_wat_thresh, NumericVector param_BASEFLOW_sur_k, NumericVector param_Lake_Eva_vic_gamma, NumericVector param_Lake_acp_storeFactor, NumericVector param_Lake_acp_gamma, NumericVector param_Riverlak_Eva_vic_gamma, NumericVector param_Riverlak_lin_storeFactor, bool if_allVariExport = false) {
        typedef SEXP(*Ptr_WaterGAP3_HL)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_WaterGAP3_HL p_WaterGAP3_HL = NULL;
        if (p_WaterGAP3_HL == NULL) {
            validateSignature("List(*WaterGAP3_HL)(int,int,NumericMatrix,NumericMatrix,NumericMatrix,IntegerVector,NumericMatrix,NumericVector,NumericVector,NumericVector,NumericMatrix,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,List,List,IntegerVector,NumericVector,NumericVector,NumericVector,IntegerVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,LogicalVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,bool)");
            p_WaterGAP3_HL = (Ptr_WaterGAP3_HL)R_GetCCallable("WaterGAP3", "_WaterGAP3_WaterGAP3_HL");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_WaterGAP3_HL(Shield<SEXP>(Rcpp::wrap(n_time)), Shield<SEXP>(Rcpp::wrap(n_spat)), Shield<SEXP>(Rcpp::wrap(ATMOS_precipitation_mm)), Shield<SEXP>(Rcpp::wrap(ATMOS_temperature_Cel)), Shield<SEXP>(Rcpp::wrap(ATMOS_potentialEvatrans_mm)), Shield<SEXP>(Rcpp::wrap(Upstream_cellNumber_int)), Shield<SEXP>(Rcpp::wrap(Upstream_streamflow_m3)), Shield<SEXP>(Rcpp::wrap(SNOW_ice_mm)), Shield<SEXP>(Rcpp::wrap(LAND_builtRatio_1)), Shield<SEXP>(Rcpp::wrap(LAND_interceptWater_mm)), Shield<SEXP>(Rcpp::wrap(LAND_interceptCapacity_mm)), Shield<SEXP>(Rcpp::wrap(SOIL_water_mm)), Shield<SEXP>(Rcpp::wrap(SOIL_capacity_mm)), Shield<SEXP>(Rcpp::wrap(SOIL_potentialPercola_mm)), Shield<SEXP>(Rcpp::wrap(GROUND_water_mm)), Shield<SEXP>(Rcpp::wrap(RIVER_water_m3)), Shield<SEXP>(Rcpp::wrap(RIVER_length_km)), Shield<SEXP>(Rcpp::wrap(RIVER_velocity_km)), Shield<SEXP>(Rcpp::wrap(CELL_landArea_km2)), Shield<SEXP>(Rcpp::wrap(CELL_cellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(CELL_inflowCellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(Lake_cellNumber_int)), Shield<SEXP>(Rcpp::wrap(Lake_water_m3)), Shield<SEXP>(Rcpp::wrap(Lake_area_km2)), Shield<SEXP>(Rcpp::wrap(Lake_capacity_m3)), Shield<SEXP>(Rcpp::wrap(Riverlak_cellNumber_int)), Shield<SEXP>(Rcpp::wrap(Riverlak_water_m3)), Shield<SEXP>(Rcpp::wrap(Riverlak_area_km2)), Shield<SEXP>(Rcpp::wrap(Riverlak_capacity_m3)), Shield<SEXP>(Rcpp::wrap(param_ATMOS_thr_Ts)), Shield<SEXP>(Rcpp::wrap(param_SNOW_fac_f)), Shield<SEXP>(Rcpp::wrap(param_SNOW_fac_Tmelt)), Shield<SEXP>(Rcpp::wrap(param_EVATRANS_sup_k)), Shield<SEXP>(Rcpp::wrap(param_EVATRANS_sup_gamma)), Shield<SEXP>(Rcpp::wrap(param_EVATRANS_wat_petmax)), Shield<SEXP>(Rcpp::wrap(param_INFILT_hbv_beta)), Shield<SEXP>(Rcpp::wrap(param_PERCOLA_wat_01)), Shield<SEXP>(Rcpp::wrap(param_PERCOLA_wat_k)), Shield<SEXP>(Rcpp::wrap(param_PERCOLA_wat_thresh)), Shield<SEXP>(Rcpp::wrap(param_BASEFLOW_sur_k)), Shield<SEXP>(Rcpp::wrap(param_Lake_Eva_vic_gamma)), Shield<SEXP>(Rcpp::wrap(param_Lake_acp_storeFactor)), Shield<SEXP>(Rcpp::wrap(param_Lake_acp_gamma)), Shield<SEXP>(Rcpp::wrap(param_Riverlak_Eva_vic_gamma)), Shield<SEXP>(Rcpp::wrap(param_Riverlak_lin_storeFactor)), Shield<SEXP>(Rcpp::wrap(if_allVariExport)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline void write_nc_WG3(NumericMatrix mat_Data_WG3, IntegerVector dim_Time, IntegerVector dim_Spat, std::string path_File, std::string name_Variable, std::string str_Continent, std::string suffix_File) {
        typedef SEXP(*Ptr_write_nc_WG3)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_write_nc_WG3 p_write_nc_WG3 = NULL;
        if (p_write_nc_WG3 == NULL) {
            validateSignature("void(*write_nc_WG3)(NumericMatrix,IntegerVector,IntegerVector,std::string,std::string,std::string,std::string)");
            p_write_nc_WG3 = (Ptr_write_nc_WG3)R_GetCCallable("WaterGAP3", "_WaterGAP3_write_nc_WG3");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_write_nc_WG3(Shield<SEXP>(Rcpp::wrap(mat_Data_WG3)), Shield<SEXP>(Rcpp::wrap(dim_Time)), Shield<SEXP>(Rcpp::wrap(dim_Spat)), Shield<SEXP>(Rcpp::wrap(path_File)), Shield<SEXP>(Rcpp::wrap(name_Variable)), Shield<SEXP>(Rcpp::wrap(str_Continent)), Shield<SEXP>(Rcpp::wrap(suffix_File)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
    }

    inline NumericMatrix read_nc_WG3(std::string path_File, std::string name_Variable) {
        typedef SEXP(*Ptr_read_nc_WG3)(SEXP,SEXP);
        static Ptr_read_nc_WG3 p_read_nc_WG3 = NULL;
        if (p_read_nc_WG3 == NULL) {
            validateSignature("NumericMatrix(*read_nc_WG3)(std::string,std::string)");
            p_read_nc_WG3 = (Ptr_read_nc_WG3)R_GetCCallable("WaterGAP3", "_WaterGAP3_read_nc_WG3");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_read_nc_WG3(Shield<SEXP>(Rcpp::wrap(path_File)), Shield<SEXP>(Rcpp::wrap(name_Variable)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

    inline Rcpp::IntegerVector read_nc_dim_WG3(std::string path_File, std::string dim_name) {
        typedef SEXP(*Ptr_read_nc_dim_WG3)(SEXP,SEXP);
        static Ptr_read_nc_dim_WG3 p_read_nc_dim_WG3 = NULL;
        if (p_read_nc_dim_WG3 == NULL) {
            validateSignature("Rcpp::IntegerVector(*read_nc_dim_WG3)(std::string,std::string)");
            p_read_nc_dim_WG3 = (Ptr_read_nc_dim_WG3)R_GetCCallable("WaterGAP3", "_WaterGAP3_read_nc_dim_WG3");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_read_nc_dim_WG3(Shield<SEXP>(Rcpp::wrap(path_File)), Shield<SEXP>(Rcpp::wrap(dim_name)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::IntegerVector >(rcpp_result_gen);
    }

}

#endif // RCPP_WaterGAP3_RCPPEXPORTS_H_GEN_
