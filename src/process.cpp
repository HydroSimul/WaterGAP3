#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' Hydrological Process
//' @name process
//' @inheritParams all_vari
//' @return atmos_snow_mm (mm/m2/TS) snowfall volume
//' @param param_atmos_thr_Ts <-1, 3> (Cel) threshold air temperature that snow, parameter for [atmosSnow_ThresholdT()]
//' @export
// [[Rcpp::export]]
NumericVector atmosSnow_ThresholdT(
   NumericVector atmos_precipitation_mm,
   NumericVector atmos_temperature_Cel,
   NumericVector param_atmos_thr_Ts
)
{
 return ifelse(atmos_temperature_Cel > param_atmos_thr_Ts, 0, atmos_precipitation_mm);
}


//' @rdname process
//' @param param_evatrans_tur_k <0.6, 1> parameter for [evatransPotential_TurcWendling()], higher value when closer to the sea
//' @return potential evapotranspiration (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector evatransPotential_TurcWendling(
   NumericVector atmos_temperature_Cel,
   NumericVector atmos_solarRadiat_MJ,
   NumericVector param_evatrans_tur_k
)
{
 return pmax((atmos_solarRadiat_MJ * 100 + 3.875 * 24 * param_evatrans_tur_k) * (atmos_temperature_Cel + 22) / 150 / (atmos_temperature_Cel + 123), 0);
 // return (atmos_solarRadiat_MJ * 100 + 3.875 * time_step_h * param_evatrans_tur_k) * (atmos_temperature_Cel + 22) / 150 / (atmos_temperature_Cel + 123);
}


//' @rdname process
//' @param param_evatrans_vic_gamma <0.2, 5> parameter for [evatransActual_VIC()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_VIC(
   NumericVector atmos_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_evatrans_vic_gamma
)
{
 NumericVector AET, k_;

 k_ = 1 - vecpow((1- water_mm / capacity_mm), param_evatrans_vic_gamma);
 AET = atmos_potentialEvatrans_mm * k_;
 return ifelse(AET > water_mm, water_mm, AET);
}


//' @rdname process
//' @param param_evatrans_ubc_gamma <0.5, 2> parameter for [evatransActual_UBC()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_UBC(
   NumericVector atmos_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_evatrans_ubc_gamma
)
{
 NumericVector diff_mm, AET, k_;
 diff_mm = capacity_mm - water_mm;



 k_ = vecpow10(- diff_mm / (param_evatrans_ubc_gamma * capacity_mm));
 AET = atmos_potentialEvatrans_mm * k_;
 return ifelse(AET > water_mm, water_mm, AET);
}

//' @rdname process
//' @param param_snow_fac_Tmelt <0, 3> (Cel) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_snow_fac_f <0.05, 2> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
NumericVector snowMelt_Factor(
   NumericVector snow_ice_mm,
   NumericVector atmos_temperature_Cel,
   NumericVector param_snow_fac_f,
   NumericVector param_snow_fac_Tmelt
)
{
 NumericVector diff_T, snow_melt_mm;
 diff_T = atmos_temperature_Cel - param_snow_fac_Tmelt;
 diff_T = ifelse(diff_T > 0, diff_T, 0);

 snow_melt_mm = param_snow_fac_f * 24 * diff_T;
 // snow_melt_mm = param_snow_fac_f * time_step_h * diff_T;
 return ifelse(snow_melt_mm > snow_ice_mm, snow_ice_mm, snow_melt_mm) ;
}

//' @rdname process
//' @param param_infilt_hbv_beta <0.001, 5> parameters for [infilt_HBV()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_HBV(
   NumericVector land_water_mm,
   NumericVector soil_water_mm,
   NumericVector soil_capacity_mm,
   NumericVector param_infilt_hbv_beta
)
{
 NumericVector soil_diff_mm, infilt_water_mm, k_, limit_mm;

 soil_diff_mm = soil_capacity_mm - soil_water_mm;
 limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);

 k_ = (1 - vecpow(soil_water_mm / soil_capacity_mm, param_infilt_hbv_beta));

 infilt_water_mm = land_water_mm * k_;

 return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}


//' @rdname process
//' @param param_percola_arn_thresh <0.1, 0.9> coefficient parameter for [percola_ThreshPow()]
//' @param param_percola_arn_k <0.1, 1> exponential parameter for [percola_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector percola_Arno(
   NumericVector soil_water_mm,
   NumericVector soil_capacity_mm,
   NumericVector soil_potentialPercola_mm,
   NumericVector param_percola_arn_thresh,
   NumericVector param_percola_arn_k
)
{
 NumericVector percola_, percola_1, percola_2, Ws_Wc;
 Ws_Wc = soil_capacity_mm * param_percola_arn_thresh;
 percola_1 = param_percola_arn_k * soil_potentialPercola_mm / (soil_capacity_mm) * soil_water_mm;
 percola_2 = param_percola_arn_k * soil_potentialPercola_mm / (soil_capacity_mm) * soil_water_mm + soil_potentialPercola_mm * (1 - param_percola_arn_k) * pow((soil_water_mm - Ws_Wc) / (soil_capacity_mm - Ws_Wc),2);
 percola_ = ifelse(soil_water_mm < Ws_Wc, percola_1, percola_2);
 percola_ = ifelse(soil_potentialPercola_mm > Ws_Wc, soil_water_mm, percola_);
 percola_ = ifelse(percola_ > soil_potentialPercola_mm, soil_potentialPercola_mm, percola_);
 return ifelse(percola_ > soil_water_mm, soil_water_mm, percola_) ;
}

//' @rdname process
//' @param param_baseflow_grf_gamma <2, 7> exponential parameter for [baseflow_GR4Jfix()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_GR4Jfix(
   NumericVector ground_water_mm,
   NumericVector ground_capacity_mm,
   NumericVector param_baseflow_grf_gamma
)
{
 NumericVector baseflow_, k_;

 k_ = 1 - vecpow((1 + vecpow(ground_water_mm / ground_capacity_mm, param_baseflow_grf_gamma)), -1.0 / param_baseflow_grf_gamma);
 baseflow_ = k_ * ground_water_mm;

 return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}

//' @rdname process
//' @param param_lake_acp_storeFactor <uknow> parameter for [lake_AcceptPow()],
//' @param param_lake_acp_gamma <uknow> parameter for [lake_AcceptPow()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector lake_AcceptPow(
   NumericVector lake_water_m3,
   NumericVector lake_inflow_m3,
   NumericVector lake_capacity_m3,
   NumericVector param_lake_acp_storeFactor,
   NumericVector param_lake_acp_gamma
)
{
 lake_water_m3 += lake_inflow_m3;

 NumericVector lake_outflow_m3 = (1 / param_lake_acp_storeFactor) * vecpow(pmin(lake_water_m3 / lake_capacity_m3, 1), param_lake_acp_gamma);
 lake_water_m3 += -lake_outflow_m3;

 NumericVector lake_overflow_m3 = pmax(lake_water_m3 -  lake_capacity_m3, 0);
 lake_water_m3 = pmin(lake_water_m3, lake_capacity_m3);

 lake_outflow_m3 += lake_overflow_m3;

 return (lake_outflow_m3);
}


//' @rdname process
//' @export
// [[Rcpp::export]]
NumericVector river_LinearResorvoir(
   NumericVector river_water_m3,
   NumericVector river_inflow_m3,
   NumericVector river_velocity_km,
   NumericVector river_length_km
)
{

 NumericVector river_paramK_TS = pmax(river_length_km / river_velocity_km, 1.);

 // return river_water_m3 * (1 / (river_paramK_TS + 0.5)) + river_inflow_m3 * (0.5 / (river_paramK_TS + 0.5));
 return river_water_m3 * (1 - exp(-1. / river_paramK_TS)) + river_inflow_m3 * (1 - river_paramK_TS * (1 - exp(-1. / river_paramK_TS)));
}

//' @rdname process
//' @param param_riverlake_lin_storeFactor <uknow> parameter for [riverlake_LinearResorvoir()],
//' @export
// [[Rcpp::export]]
NumericVector riverlake_LinearResorvoir(
   NumericVector riverlake_water_m3,
   NumericVector riverlake_inflow_m3,
   NumericVector riverlake_capacity_m3,
   NumericVector param_riverlake_lin_storeFactor
)
{
 riverlake_water_m3 += riverlake_inflow_m3;

 NumericVector riverlake_outflow_m3 = riverlake_water_m3 * (1 - exp(-1. / param_riverlake_lin_storeFactor)) + riverlake_inflow_m3 * (1 - param_riverlake_lin_storeFactor * (1 - exp(-1. / param_riverlake_lin_storeFactor)));

 NumericVector riverlake_overflow_m3 = pmax(riverlake_water_m3 -  riverlake_capacity_m3, 0);
 riverlake_water_m3 = pmin(riverlake_water_m3, riverlake_capacity_m3);

 riverlake_outflow_m3 += riverlake_overflow_m3;

 return (riverlake_outflow_m3);
}

//' @rdname process
//' @param param_reservoir_han_alpha <uknow> parameter for [reservoir_Hanasaki()],
//' @param param_reservoir_han_kDemand <uknow> parameter for [reservoir_Hanasaki()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector reservoir_Hanasaki(
   NumericVector reservoir_water_m3,
   NumericVector reservoir_inflow_m3,
   NumericVector reservoir_capacity_m3,
   NumericVector reservoir_demand_m3,
   NumericVector reservoir_yearInflow_m3,
   NumericVector reservoir_yearDemand_m3,
   NumericVector &reservoir_yearRelase_m3,
   LogicalVector reservoir_isOperateStart_01,
   LogicalVector reservoir_isIrrigate_01,
   NumericVector param_reservoir_han_alpha,
   NumericVector param_reservoir_han_kDemand
)
{
 reservoir_water_m3 += reservoir_inflow_m3;

 NumericVector reservoir_releaseCoefficient_ = reservoir_water_m3 / (param_reservoir_han_alpha * reservoir_capacity_m3); // eq-3
 reservoir_releaseCoefficient_ = pmax(reservoir_releaseCoefficient_, 0.1);
 reservoir_yearRelase_m3 = ifelse(reservoir_isOperateStart_01, reservoir_releaseCoefficient_ * reservoir_yearInflow_m3, reservoir_yearRelase_m3); // eq-2

 NumericVector reservoir_releaseProvis_m3 = reservoir_yearInflow_m3; // eq-4
 reservoir_releaseProvis_m3 = ifelse(reservoir_isIrrigate_01,
                                     0.5 * reservoir_yearInflow_m3 * (1 + param_reservoir_han_kDemand * reservoir_demand_m3 / reservoir_yearDemand_m3),
                                     reservoir_releaseProvis_m3);  // eq-5
 reservoir_releaseProvis_m3 = ifelse(reservoir_yearDemand_m3 < 0.5 * reservoir_yearInflow_m3,
                                     reservoir_yearInflow_m3 + param_reservoir_han_kDemand * reservoir_demand_m3 - reservoir_yearDemand_m3,
                                     reservoir_releaseProvis_m3);  // eq-5
 NumericVector reservoir_inflowRatio_ = reservoir_capacity_m3 / reservoir_yearInflow_m3; // eq-7
 NumericVector temp_inflowRatio_ = 4 * reservoir_inflowRatio_ * reservoir_inflowRatio_; // eq-7

 return ifelse(reservoir_inflowRatio_ > 0.5,
               reservoir_releaseCoefficient_ * reservoir_releaseProvis_m3,
               temp_inflowRatio_ * reservoir_releaseCoefficient_ * reservoir_releaseProvis_m3 + (1 - temp_inflowRatio_) * reservoir_inflow_m3); //  eq-7
}


NumericVector inflow_add(
    NumericVector num_Outflow_LastStep,
    IntegerMatrix int_InflowCell
)
{

  int n_Cell = int_InflowCell.nrow();  // Number of rows (cells)
  NumericVector num_Inflow_m3(n_Cell, 0.); // Initialize the result vector with NA

  for (int i_Cell = 0; i_Cell < n_Cell; i_Cell++) {
    double inflow_sum = 0.0;

    for (int i_InflowCell = 0; i_InflowCell < int_InflowCell.ncol(); i_InflowCell++) {
      int index = int_InflowCell(i_Cell, i_InflowCell);
      if (index == NA_INTEGER) break;
      inflow_sum += num_Outflow_LastStep[index - 1];
    }
    num_Inflow_m3[i_Cell] = inflow_sum;
  }

  return num_Inflow_m3;
}



IntegerVector get_idx_step(IntegerVector int_Cell, IntegerVector int_Step) {
  // Find the match between int_Step and int_Cell
  IntegerVector int_Match = match(int_Cell, int_Step);

  return na_omit(int_Match);
}

IntegerVector get_idx_cell(IntegerVector int_Cell, IntegerVector int_Step) {
  // Find the match between int_Step and int_Cell
  IntegerVector int_Match = match(int_Step, int_Cell);

  return na_omit(int_Match);
}


//' @rdname process
//' @export
// [[Rcpp::export]]
NumericVector confluen_WaterGAP3(
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
   NumericVector param_lake_acp_storeFactor,
   NumericVector param_lake_acp_gamma,
   NumericVector param_riverlake_lin_storeFactor,
   NumericVector param_reservoir_han_alpha,
   NumericVector param_reservoir_han_kDemand
)
{

 int n_Cell = confluen_cellInflow_m3.size();
 NumericVector confluen_outflow_m3(n_Cell), waterbody_Evatrans_mm, step_RiverOutflow_m3,
 step_RiverlakeOutflow_m3, step_ReservoirOutflow_m3,
 lake_Outflow_m3, lake_Evatrans_mm,
 step_Inflow_Riverlake_m3;

 IntegerVector idx_Cell_Step,
 idx_Riverlake_Step, idx_Step_Riverlake,
 idx_Reservoir_Step, idx_Step_Reservoir;
 int n_Step = basin_cellNumberStep_int.size();

 // Step i later with Inflow
 for (int i_Step = 0; i_Step < n_Step; i_Step++)
 {

   idx_Cell_Step = basin_cellNumberStep_int[i_Step];
   NumericVector step_UpstreamInflow_m3(idx_Cell_Step.size(), 0.);
   // Inflow upstream
   if (i_Step > 0) {

     step_UpstreamInflow_m3 = inflow_add(
       confluen_outflow_m3,
       basin_inflowCellNumberStep_int[i_Step]
     );

   }


   // river segment
   NumericVector step_CellInflow = subset_get(confluen_cellInflow_m3, idx_Cell_Step),
     step_RiverWater= subset_get(river_water_m3, idx_Cell_Step),
     step_RiverInflow = step_UpstreamInflow_m3 + step_CellInflow;

   step_RiverOutflow_m3 = river_LinearResorvoir(
     step_RiverWater,
     step_RiverInflow,
     subset_get(river_velocity_km, idx_Cell_Step),
     subset_get(river_length_km, idx_Cell_Step)
   );
   NumericVector step_RiverInOut = pmax(step_RiverInflow - step_RiverOutflow_m3, 0.0);
   subset_put(confluen_outflow_m3, idx_Cell_Step, step_RiverOutflow_m3);
   subset_put(river_water_m3, idx_Cell_Step, step_RiverWater + step_RiverInOut);

   // global lake (riverlake)
   idx_Riverlake_Step = get_idx_cell(riverlake_cellNumber_int, idx_Cell_Step);
   idx_Step_Riverlake = get_idx_step(riverlake_cellNumber_int, idx_Cell_Step);
   if (idx_Riverlake_Step.size() > 0) {
     NumericVector step_RiverlakeWater= subset_get(riverlake_water_m3, idx_Riverlake_Step),
       step_RiverlakeInflow = subset_get(step_UpstreamInflow_m3 + step_CellInflow, idx_Step_Riverlake);


     step_RiverlakeOutflow_m3 = riverlake_LinearResorvoir(
       step_RiverlakeWater,
       step_RiverlakeInflow,
       subset_get(riverlake_capacity_m3, idx_Riverlake_Step),
       subset_get(param_riverlake_lin_storeFactor, idx_Riverlake_Step)
     );
     NumericVector step_RiverlakeInOut = pmax(step_RiverlakeInflow - step_RiverlakeOutflow_m3, 0.0);
     subset_put(confluen_outflow_m3, idx_Riverlake_Step, step_RiverlakeOutflow_m3);
     subset_put(riverlake_water_m3, idx_Riverlake_Step, step_RiverlakeWater + step_RiverlakeInOut);
   }


   // Reservior
   idx_Reservoir_Step = get_idx_cell(reservoir_cellNumber_int, idx_Cell_Step);
   idx_Step_Reservoir = get_idx_step(reservoir_cellNumber_int, idx_Cell_Step);
   if (idx_Reservoir_Step.size() > 0 ) {
     NumericVector step_ReservoirYearRealse = subset_get(reservoir_yearRelase_m3, idx_Reservoir_Step),
       step_ReservoirWater= subset_get(reservoir_water_m3, idx_Reservoir_Step),
       step_ReservoirInflow = subset_get(step_UpstreamInflow_m3 + step_CellInflow, idx_Step_Reservoir);
     step_ReservoirOutflow_m3 = reservoir_Hanasaki(
       step_ReservoirWater,
       step_ReservoirInflow,
       subset_get(reservoir_capacity_m3, idx_Reservoir_Step),
       subset_get(reservoir_demand_m3, idx_Reservoir_Step),
       subset_get(reservoir_yearInflow_m3, idx_Reservoir_Step),
       subset_get(reservoir_yearDemand_m3, idx_Reservoir_Step),
       step_ReservoirYearRealse,
       subset_get_logical(reservoir_isOperateStart_01, idx_Reservoir_Step),
       subset_get_logical(reservoir_isIrrigate_01, idx_Reservoir_Step),
       subset_get(param_reservoir_han_alpha, idx_Reservoir_Step),
       subset_get(param_reservoir_han_kDemand, idx_Reservoir_Step)
     );
     NumericVector step_ReservoirInOut = pmax(step_ReservoirInflow - step_ReservoirOutflow_m3, 0.0);
     subset_put(reservoir_yearRelase_m3, idx_Reservoir_Step, step_ReservoirYearRealse);
     subset_put(confluen_outflow_m3, idx_Reservoir_Step, step_ReservoirOutflow_m3);
     subset_put(reservoir_water_m3, idx_Reservoir_Step, step_ReservoirWater + step_ReservoirInOut);
   }


 }


 return confluen_outflow_m3;

}





