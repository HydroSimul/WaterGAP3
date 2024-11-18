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

    inline NumericVector subset_get(NumericVector vec_Data, IntegerVector int_Index) {
        typedef SEXP(*Ptr_subset_get)(SEXP,SEXP);
        static Ptr_subset_get p_subset_get = NULL;
        if (p_subset_get == NULL) {
            validateSignature("NumericVector(*subset_get)(NumericVector,IntegerVector)");
            p_subset_get = (Ptr_subset_get)R_GetCCallable("WaterGAP3", "_WaterGAP3_subset_get");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_subset_get(Shield<SEXP>(Rcpp::wrap(vec_Data)), Shield<SEXP>(Rcpp::wrap(int_Index)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline LogicalVector subset_get_logical(LogicalVector vec_Data, IntegerVector int_Index) {
        typedef SEXP(*Ptr_subset_get_logical)(SEXP,SEXP);
        static Ptr_subset_get_logical p_subset_get_logical = NULL;
        if (p_subset_get_logical == NULL) {
            validateSignature("LogicalVector(*subset_get_logical)(LogicalVector,IntegerVector)");
            p_subset_get_logical = (Ptr_subset_get_logical)R_GetCCallable("WaterGAP3", "_WaterGAP3_subset_get_logical");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_subset_get_logical(Shield<SEXP>(Rcpp::wrap(vec_Data)), Shield<SEXP>(Rcpp::wrap(int_Index)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<LogicalVector >(rcpp_result_gen);
    }

    inline void subset_put(NumericVector& vec_Data, IntegerVector int_Index, NumericVector vec_DataPut) {
        typedef SEXP(*Ptr_subset_put)(SEXP,SEXP,SEXP);
        static Ptr_subset_put p_subset_put = NULL;
        if (p_subset_put == NULL) {
            validateSignature("void(*subset_put)(NumericVector&,IntegerVector,NumericVector)");
            p_subset_put = (Ptr_subset_put)R_GetCCallable("WaterGAP3", "_WaterGAP3_subset_put");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_subset_put(Shield<SEXP>(Rcpp::wrap(vec_Data)), Shield<SEXP>(Rcpp::wrap(int_Index)), Shield<SEXP>(Rcpp::wrap(vec_DataPut)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
    }

    inline NumericVector atmosSnow_ThresholdT(NumericVector AtmoS_precipitation_mm, NumericVector AtmoS_temperature_Cel, NumericVector param_atmos_thr_Ts) {
        typedef SEXP(*Ptr_atmosSnow_ThresholdT)(SEXP,SEXP,SEXP);
        static Ptr_atmosSnow_ThresholdT p_atmosSnow_ThresholdT = NULL;
        if (p_atmosSnow_ThresholdT == NULL) {
            validateSignature("NumericVector(*atmosSnow_ThresholdT)(NumericVector,NumericVector,NumericVector)");
            p_atmosSnow_ThresholdT = (Ptr_atmosSnow_ThresholdT)R_GetCCallable("WaterGAP3", "_WaterGAP3_atmosSnow_ThresholdT");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_atmosSnow_ThresholdT(Shield<SEXP>(Rcpp::wrap(AtmoS_precipitation_mm)), Shield<SEXP>(Rcpp::wrap(AtmoS_temperature_Cel)), Shield<SEXP>(Rcpp::wrap(param_atmos_thr_Ts)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector evatransPotential_TurcWendling(NumericVector AtmoS_temperature_Cel, NumericVector AtmoS_solarRadiat_MJ, NumericVector param_evatrans_tur_k) {
        typedef SEXP(*Ptr_evatransPotential_TurcWendling)(SEXP,SEXP,SEXP);
        static Ptr_evatransPotential_TurcWendling p_evatransPotential_TurcWendling = NULL;
        if (p_evatransPotential_TurcWendling == NULL) {
            validateSignature("NumericVector(*evatransPotential_TurcWendling)(NumericVector,NumericVector,NumericVector)");
            p_evatransPotential_TurcWendling = (Ptr_evatransPotential_TurcWendling)R_GetCCallable("WaterGAP3", "_WaterGAP3_evatransPotential_TurcWendling");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_evatransPotential_TurcWendling(Shield<SEXP>(Rcpp::wrap(AtmoS_temperature_Cel)), Shield<SEXP>(Rcpp::wrap(AtmoS_solarRadiat_MJ)), Shield<SEXP>(Rcpp::wrap(param_evatrans_tur_k)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector evatransActual_VIC(NumericVector AtmoS_potentialEvatrans_mm, NumericVector water_mm, NumericVector capacity_mm, NumericVector param_evatrans_vic_gamma) {
        typedef SEXP(*Ptr_evatransActual_VIC)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_evatransActual_VIC p_evatransActual_VIC = NULL;
        if (p_evatransActual_VIC == NULL) {
            validateSignature("NumericVector(*evatransActual_VIC)(NumericVector,NumericVector,NumericVector,NumericVector)");
            p_evatransActual_VIC = (Ptr_evatransActual_VIC)R_GetCCallable("WaterGAP3", "_WaterGAP3_evatransActual_VIC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_evatransActual_VIC(Shield<SEXP>(Rcpp::wrap(AtmoS_potentialEvatrans_mm)), Shield<SEXP>(Rcpp::wrap(water_mm)), Shield<SEXP>(Rcpp::wrap(capacity_mm)), Shield<SEXP>(Rcpp::wrap(param_evatrans_vic_gamma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector evatransActual_UBC(NumericVector AtmoS_potentialEvatrans_mm, NumericVector water_mm, NumericVector capacity_mm, NumericVector param_evatrans_ubc_gamma) {
        typedef SEXP(*Ptr_evatransActual_UBC)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_evatransActual_UBC p_evatransActual_UBC = NULL;
        if (p_evatransActual_UBC == NULL) {
            validateSignature("NumericVector(*evatransActual_UBC)(NumericVector,NumericVector,NumericVector,NumericVector)");
            p_evatransActual_UBC = (Ptr_evatransActual_UBC)R_GetCCallable("WaterGAP3", "_WaterGAP3_evatransActual_UBC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_evatransActual_UBC(Shield<SEXP>(Rcpp::wrap(AtmoS_potentialEvatrans_mm)), Shield<SEXP>(Rcpp::wrap(water_mm)), Shield<SEXP>(Rcpp::wrap(capacity_mm)), Shield<SEXP>(Rcpp::wrap(param_evatrans_ubc_gamma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector snowMelt_Factor(NumericVector snoW_ice_mm, NumericVector AtmoS_temperature_Cel, NumericVector param_snow_fac_f, NumericVector param_snow_fac_Tmelt) {
        typedef SEXP(*Ptr_snowMelt_Factor)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_snowMelt_Factor p_snowMelt_Factor = NULL;
        if (p_snowMelt_Factor == NULL) {
            validateSignature("NumericVector(*snowMelt_Factor)(NumericVector,NumericVector,NumericVector,NumericVector)");
            p_snowMelt_Factor = (Ptr_snowMelt_Factor)R_GetCCallable("WaterGAP3", "_WaterGAP3_snowMelt_Factor");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_snowMelt_Factor(Shield<SEXP>(Rcpp::wrap(snoW_ice_mm)), Shield<SEXP>(Rcpp::wrap(AtmoS_temperature_Cel)), Shield<SEXP>(Rcpp::wrap(param_snow_fac_f)), Shield<SEXP>(Rcpp::wrap(param_snow_fac_Tmelt)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector infilt_HBV(NumericVector lanD_water_mm, NumericVector soiL_water_mm, NumericVector soiL_capacity_mm, NumericVector param_infilt_hbv_beta) {
        typedef SEXP(*Ptr_infilt_HBV)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_infilt_HBV p_infilt_HBV = NULL;
        if (p_infilt_HBV == NULL) {
            validateSignature("NumericVector(*infilt_HBV)(NumericVector,NumericVector,NumericVector,NumericVector)");
            p_infilt_HBV = (Ptr_infilt_HBV)R_GetCCallable("WaterGAP3", "_WaterGAP3_infilt_HBV");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_infilt_HBV(Shield<SEXP>(Rcpp::wrap(lanD_water_mm)), Shield<SEXP>(Rcpp::wrap(soiL_water_mm)), Shield<SEXP>(Rcpp::wrap(soiL_capacity_mm)), Shield<SEXP>(Rcpp::wrap(param_infilt_hbv_beta)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector percola_Arno(NumericVector soiL_water_mm, NumericVector soiL_capacity_mm, NumericVector soiL_potentialPercola_mm, NumericVector param_percola_arn_thresh, NumericVector param_percola_arn_k) {
        typedef SEXP(*Ptr_percola_Arno)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_percola_Arno p_percola_Arno = NULL;
        if (p_percola_Arno == NULL) {
            validateSignature("NumericVector(*percola_Arno)(NumericVector,NumericVector,NumericVector,NumericVector,NumericVector)");
            p_percola_Arno = (Ptr_percola_Arno)R_GetCCallable("WaterGAP3", "_WaterGAP3_percola_Arno");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_percola_Arno(Shield<SEXP>(Rcpp::wrap(soiL_water_mm)), Shield<SEXP>(Rcpp::wrap(soiL_capacity_mm)), Shield<SEXP>(Rcpp::wrap(soiL_potentialPercola_mm)), Shield<SEXP>(Rcpp::wrap(param_percola_arn_thresh)), Shield<SEXP>(Rcpp::wrap(param_percola_arn_k)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector baseflow_GR4Jfix(NumericVector grounD_water_mm, NumericVector grounD_capacity_mm, NumericVector param_baseflow_grf_gamma) {
        typedef SEXP(*Ptr_baseflow_GR4Jfix)(SEXP,SEXP,SEXP);
        static Ptr_baseflow_GR4Jfix p_baseflow_GR4Jfix = NULL;
        if (p_baseflow_GR4Jfix == NULL) {
            validateSignature("NumericVector(*baseflow_GR4Jfix)(NumericVector,NumericVector,NumericVector)");
            p_baseflow_GR4Jfix = (Ptr_baseflow_GR4Jfix)R_GetCCallable("WaterGAP3", "_WaterGAP3_baseflow_GR4Jfix");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_baseflow_GR4Jfix(Shield<SEXP>(Rcpp::wrap(grounD_water_mm)), Shield<SEXP>(Rcpp::wrap(grounD_capacity_mm)), Shield<SEXP>(Rcpp::wrap(param_baseflow_grf_gamma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector lake_AcceptPow(NumericVector lake_water_m3, NumericVector lake_inflow_m3, NumericVector lake_capacity_m3, NumericVector param_lake_acp_storeFactor, NumericVector param_lake_acp_gamma) {
        typedef SEXP(*Ptr_lake_AcceptPow)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_lake_AcceptPow p_lake_AcceptPow = NULL;
        if (p_lake_AcceptPow == NULL) {
            validateSignature("NumericVector(*lake_AcceptPow)(NumericVector,NumericVector,NumericVector,NumericVector,NumericVector)");
            p_lake_AcceptPow = (Ptr_lake_AcceptPow)R_GetCCallable("WaterGAP3", "_WaterGAP3_lake_AcceptPow");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_lake_AcceptPow(Shield<SEXP>(Rcpp::wrap(lake_water_m3)), Shield<SEXP>(Rcpp::wrap(lake_inflow_m3)), Shield<SEXP>(Rcpp::wrap(lake_capacity_m3)), Shield<SEXP>(Rcpp::wrap(param_lake_acp_storeFactor)), Shield<SEXP>(Rcpp::wrap(param_lake_acp_gamma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector river_LinearResorvoir(NumericVector riveR_water_m3, NumericVector riveR_inflow_m3, NumericVector riveR_velocity_km, NumericVector riveR_length_km) {
        typedef SEXP(*Ptr_river_LinearResorvoir)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_river_LinearResorvoir p_river_LinearResorvoir = NULL;
        if (p_river_LinearResorvoir == NULL) {
            validateSignature("NumericVector(*river_LinearResorvoir)(NumericVector,NumericVector,NumericVector,NumericVector)");
            p_river_LinearResorvoir = (Ptr_river_LinearResorvoir)R_GetCCallable("WaterGAP3", "_WaterGAP3_river_LinearResorvoir");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_river_LinearResorvoir(Shield<SEXP>(Rcpp::wrap(riveR_water_m3)), Shield<SEXP>(Rcpp::wrap(riveR_inflow_m3)), Shield<SEXP>(Rcpp::wrap(riveR_velocity_km)), Shield<SEXP>(Rcpp::wrap(riveR_length_km)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector riverlake_LinearResorvoir(NumericVector riverlake_water_m3, NumericVector riverlake_inflow_m3, NumericVector riverlake_capacity_m3, NumericVector param_riverlake_lin_storeFactor) {
        typedef SEXP(*Ptr_riverlake_LinearResorvoir)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_riverlake_LinearResorvoir p_riverlake_LinearResorvoir = NULL;
        if (p_riverlake_LinearResorvoir == NULL) {
            validateSignature("NumericVector(*riverlake_LinearResorvoir)(NumericVector,NumericVector,NumericVector,NumericVector)");
            p_riverlake_LinearResorvoir = (Ptr_riverlake_LinearResorvoir)R_GetCCallable("WaterGAP3", "_WaterGAP3_riverlake_LinearResorvoir");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_riverlake_LinearResorvoir(Shield<SEXP>(Rcpp::wrap(riverlake_water_m3)), Shield<SEXP>(Rcpp::wrap(riverlake_inflow_m3)), Shield<SEXP>(Rcpp::wrap(riverlake_capacity_m3)), Shield<SEXP>(Rcpp::wrap(param_riverlake_lin_storeFactor)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector reservoir_Hanasaki(NumericVector reservoir_water_m3, NumericVector reservoir_inflow_m3, NumericVector reservoir_capacity_m3, NumericVector reservoir_demand_m3, NumericVector reservoir_yearInflow_m3, NumericVector reservoir_yearDemand_m3, NumericVector& reservoir_yearRelase_m3, LogicalVector reservoir_isOperateStart_01, LogicalVector reservoir_isIrrigate_01, NumericVector param_reservoir_han_alpha, NumericVector param_reservoir_han_kDemand) {
        typedef SEXP(*Ptr_reservoir_Hanasaki)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_reservoir_Hanasaki p_reservoir_Hanasaki = NULL;
        if (p_reservoir_Hanasaki == NULL) {
            validateSignature("NumericVector(*reservoir_Hanasaki)(NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector&,LogicalVector,LogicalVector,NumericVector,NumericVector)");
            p_reservoir_Hanasaki = (Ptr_reservoir_Hanasaki)R_GetCCallable("WaterGAP3", "_WaterGAP3_reservoir_Hanasaki");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_reservoir_Hanasaki(Shield<SEXP>(Rcpp::wrap(reservoir_water_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_inflow_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_capacity_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_demand_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_yearInflow_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_yearDemand_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_yearRelase_m3)), Shield<SEXP>(Rcpp::wrap(reservoir_isOperateStart_01)), Shield<SEXP>(Rcpp::wrap(reservoir_isIrrigate_01)), Shield<SEXP>(Rcpp::wrap(param_reservoir_han_alpha)), Shield<SEXP>(Rcpp::wrap(param_reservoir_han_kDemand)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector confluen_WaterGAP3(NumericVector conflueN_cellInflow_m3, NumericVector& riveR_water_m3, NumericVector riveR_length_km, NumericVector riveR_velocity_km, List celL_cellNumberStep_int, List celL_inflowCellNumberStep_int) {
        typedef SEXP(*Ptr_confluen_WaterGAP3)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_confluen_WaterGAP3 p_confluen_WaterGAP3 = NULL;
        if (p_confluen_WaterGAP3 == NULL) {
            validateSignature("NumericVector(*confluen_WaterGAP3)(NumericVector,NumericVector&,NumericVector,NumericVector,List,List)");
            p_confluen_WaterGAP3 = (Ptr_confluen_WaterGAP3)R_GetCCallable("WaterGAP3", "_WaterGAP3_confluen_WaterGAP3");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_confluen_WaterGAP3(Shield<SEXP>(Rcpp::wrap(conflueN_cellInflow_m3)), Shield<SEXP>(Rcpp::wrap(riveR_water_m3)), Shield<SEXP>(Rcpp::wrap(riveR_length_km)), Shield<SEXP>(Rcpp::wrap(riveR_velocity_km)), Shield<SEXP>(Rcpp::wrap(celL_cellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(celL_inflowCellNumberStep_int)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector confluen_WaterGAP3_L(NumericVector conflueN_cellInflow_m3, NumericVector& riveR_water_m3, NumericVector riveR_length_km, NumericVector riveR_velocity_km, List celL_cellNumberStep_int, List celL_inflowCellNumberStep_int, IntegerVector riverlake_cellNumber_int, NumericVector& riverlake_water_m3, NumericVector riverlake_capacity_m3, NumericVector param_riverlake_lin_storeFactor) {
        typedef SEXP(*Ptr_confluen_WaterGAP3_L)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_confluen_WaterGAP3_L p_confluen_WaterGAP3_L = NULL;
        if (p_confluen_WaterGAP3_L == NULL) {
            validateSignature("NumericVector(*confluen_WaterGAP3_L)(NumericVector,NumericVector&,NumericVector,NumericVector,List,List,IntegerVector,NumericVector&,NumericVector,NumericVector)");
            p_confluen_WaterGAP3_L = (Ptr_confluen_WaterGAP3_L)R_GetCCallable("WaterGAP3", "_WaterGAP3_confluen_WaterGAP3_L");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_confluen_WaterGAP3_L(Shield<SEXP>(Rcpp::wrap(conflueN_cellInflow_m3)), Shield<SEXP>(Rcpp::wrap(riveR_water_m3)), Shield<SEXP>(Rcpp::wrap(riveR_length_km)), Shield<SEXP>(Rcpp::wrap(riveR_velocity_km)), Shield<SEXP>(Rcpp::wrap(celL_cellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(celL_inflowCellNumberStep_int)), Shield<SEXP>(Rcpp::wrap(riverlake_cellNumber_int)), Shield<SEXP>(Rcpp::wrap(riverlake_water_m3)), Shield<SEXP>(Rcpp::wrap(riverlake_capacity_m3)), Shield<SEXP>(Rcpp::wrap(param_riverlake_lin_storeFactor)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

}

#endif // RCPP_WaterGAP3_RCPPEXPORTS_H_GEN_
