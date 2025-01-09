#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' Hydrological Process
//' @name process
//' @inheritParams all_vari
//' @return ATMOS_snow_mm (mm/m2/TS) snowfall volume
//' @param param_ATMOS_thr_Ts <-1, 3> (Cel) threshold air temperature that snow, parameter for [atmosSnow_ThresholdT()]
//' @export
// [[Rcpp::export]]
NumericVector atmosSnow_ThresholdT(
   NumericVector ATMOS_precipitation_mm,
   NumericVector ATMOS_temperature_Cel,
   NumericVector param_ATMOS_thr_Ts
)
{
 return ifelse(ATMOS_temperature_Cel > param_ATMOS_thr_Ts, 0, ATMOS_precipitation_mm);
}


//' @rdname process
//' @param param_EVATRANS_tur_k <0.6, 1> parameter for [evatransPotential_TurcWendling()], higher value when closer to the sea
//' @return potential evapotranspiration (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector evatransPotential_TurcWendling(
   NumericVector ATMOS_temperature_Cel,
   NumericVector ATMOS_solarRadiat_MJ,
   NumericVector param_EVATRANS_tur_k
)
{
 return pmax((ATMOS_solarRadiat_MJ * 100 + 3.875 * 24 * param_EVATRANS_tur_k) * (ATMOS_temperature_Cel + 22) / 150 / (ATMOS_temperature_Cel + 123), 0);
 // return (ATMOS_solarRadiat_MJ * 100 + 3.875 * time_step_h * param_EVATRANS_tur_k) * (ATMOS_temperature_Cel + 22) / 150 / (ATMOS_temperature_Cel + 123);
}


//' @rdname process
//' @export
// [[Rcpp::export]]
NumericVector intercep_Full(
    NumericVector ATMOS_precipitation_mm,
    NumericVector LAND_interceptWater_mm,
    NumericVector LAND_interceptCapacity_mm
)
{
  NumericVector water_diff_mm = LAND_interceptCapacity_mm - LAND_interceptWater_mm;
  NumericVector intercep_mm = ifelse(water_diff_mm > ATMOS_precipitation_mm, ATMOS_precipitation_mm, water_diff_mm);
  return intercep_mm;
  // return ifelse(water_diff_mm < 0, 0, intercep_mm);
}


//' @rdname process
//' @param param_EVATRANS_sup_k <0.1, 1> parameter for [evatransActual_SupplyPow()], ratio of this method
//' @param param_EVATRANS_sup_gamma <.1, 5> parameter for [evatransActual_SupplyPow()], exponent of this method
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_SupplyPow(
   NumericVector ATMOS_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_EVATRANS_sup_k,
   NumericVector param_EVATRANS_sup_gamma
)
{
 NumericVector AET, k_;

 k_ = param_EVATRANS_sup_k * vecpow((water_mm / capacity_mm), param_EVATRANS_sup_gamma);
 AET = ATMOS_potentialEvatrans_mm * k_;
 return ifelse(AET > water_mm, water_mm, AET);
}



//' @rdname process
//' @param param_EVATRANS_sur_k <0.1, 1> parameter for [evatransActual_SupplyRatio()], ratio of potential ET, that is estimated as actually ET
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_SupplyRatio(
   NumericVector ATMOS_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_EVATRANS_sur_k
)
{
 NumericVector AET, k_;

 k_ = water_mm / capacity_mm * param_EVATRANS_sur_k;
 AET = ATMOS_potentialEvatrans_mm * k_;
 return ifelse(AET > water_mm, water_mm, AET);
}


//' @rdname process
//' @param param_EVATRANS_wat_petmax <10, 20> parameter for [evatransActual_WaterGAP()], 10 for arid area, 20 for humid area
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_WaterGAP(
   NumericVector ATMOS_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_EVATRANS_wat_petmax
)
{
 NumericVector AET, AET_Temp;

 AET_Temp = water_mm / capacity_mm * param_EVATRANS_wat_petmax;
 AET = ifelse(AET_Temp > ATMOS_potentialEvatrans_mm, ATMOS_potentialEvatrans_mm, AET_Temp);
 return ifelse(AET > water_mm, water_mm, AET);
}

//' @rdname process
//' @param param_EVATRANS_vic_gamma <0.2, 5> parameter for [evatransActual_VIC()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_VIC(
   NumericVector ATMOS_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_EVATRANS_vic_gamma
)
{
 NumericVector AET, k_;

 k_ = 1 - vecpow((1- water_mm / capacity_mm), param_EVATRANS_vic_gamma);
 AET = ATMOS_potentialEvatrans_mm * k_;
 return ifelse(AET > water_mm, water_mm, AET);
}


//' @rdname process
//' @param param_EVATRANS_ubc_gamma <0.5, 2> parameter for [evatransActual_UBC()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_UBC(
   NumericVector ATMOS_potentialEvatrans_mm,
   NumericVector water_mm,
   NumericVector capacity_mm,
   NumericVector param_EVATRANS_ubc_gamma
)
{
 NumericVector diff_mm, AET, k_;
 diff_mm = capacity_mm - water_mm;



 k_ = vecpow10(- diff_mm / (param_EVATRANS_ubc_gamma * capacity_mm));
 AET = ATMOS_potentialEvatrans_mm * k_;
 return ifelse(AET > water_mm, water_mm, AET);
}

//' @rdname process
//' @param param_SNOW_fac_Tmelt <0, 3> (Cel) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_SNOW_fac_f <0.05, 2> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
NumericVector snowMelt_Factor(
   NumericVector SNOW_ice_mm,
   NumericVector ATMOS_temperature_Cel,
   NumericVector param_SNOW_fac_f,
   NumericVector param_SNOW_fac_Tmelt
)
{
 NumericVector diff_T, snow_melt_mm;
 diff_T = ATMOS_temperature_Cel - param_SNOW_fac_Tmelt;
 diff_T = ifelse(diff_T > 0, diff_T, 0);

 snow_melt_mm = param_SNOW_fac_f * 24 * diff_T;
 // snow_melt_mm = param_SNOW_fac_f * time_step_h * diff_T;
 return ifelse(snow_melt_mm > SNOW_ice_mm, SNOW_ice_mm, snow_melt_mm) ;
}


//' @rdname process
//' @export
// [[Rcpp::export]]
NumericMatrix landLeafAreaIndex_WaterGAP3(NumericMatrix ATMOS_temperature_Cel,
                                                NumericMatrix ATMOS_precipitation_mm,
                                                NumericVector CELL_latitude_deg,
                                                IntegerVector LAND_growUpDay_d,
                                                NumericVector LAND_leafAreaIndexMin_,
                                                NumericVector LAND_leafAreaIndexMax_,
                                                IntegerVector Time_dayOfYear_d) {
  int n_Days = ATMOS_temperature_Cel.nrow();
  int n_Grids = ATMOS_temperature_Cel.ncol();
  NumericMatrix LAND_leafAreaIndex_(n_Days, n_Grids);
  NumericVector LAND_leafAreaRatio_(n_Days);
  NumericVector range_LAI = LAND_leafAreaIndexMax_ - LAND_leafAreaIndexMin_;
  NumericVector num_Growup(366, 1.0);
  NumericVector num_Droop(366, 0.0);

  // Initialize growth and droop vectors
  for (int i = 0; i < 30; i++) {
    num_Growup[i] = (i + 1) / 30.0;
    num_Droop[i] = (30 - i) / 30.0;
  }

  for (int g = 0; g < n_Grids; ++g) {
    IntegerVector tag_Day;
    if (CELL_latitude_deg[g] >= 0) {
      tag_Day = find_locations(Time_dayOfYear_d, 1);
    } else {
      tag_Day = find_locations(Time_dayOfYear_d, 183);
    }

    int n_Year = tag_Day.size();
    tag_Day.push_back(n_Days + 1); // Add a placeholder for end of data

    for (int y = 0; y < n_Year; ++y) {
      int start_day = tag_Day[y];
      int next_year_start = tag_Day[y + 1];

      double cumsum_Perc_Temp = 0.0;
      for (int d = start_day; d < next_year_start && d < n_Days; ++d) {
        cumsum_Perc_Temp += ATMOS_precipitation_mm(d, g);

        // Temperature check over LAND_growUpDay_d[g] days
        int start_idx = std::max(0, d - LAND_growUpDay_d[g]);
        double cum_Temp_Temp = min(NumericVector(ATMOS_temperature_Cel(_, g).begin() + start_idx,
                                                 ATMOS_temperature_Cel(_, g).begin() + d + 1));

        if (cumsum_Perc_Temp > 40 && cum_Temp_Temp > 8) {
          for (int i = d; i < next_year_start && i < n_Days; ++i) {
            LAND_leafAreaRatio_(i) = num_Growup[std::min(i - d, 365)];
          }
          break;
        }
      }

      for (int d = start_day + 182; d < next_year_start && d < n_Days; ++d) {
        // Temperature check for drooping
        int start_idx = std::max(0, d - LAND_growUpDay_d[g]);
        double cum_Temp_Temp = max(NumericVector(ATMOS_temperature_Cel(_, g).begin() + start_idx,
                                                 ATMOS_temperature_Cel(_, g).begin() + d + 1));

        if (cum_Temp_Temp < 8) {
          for (int i = d; i < next_year_start && i < n_Days; ++i) {
            LAND_leafAreaRatio_(i) = num_Droop[std::min(i - d, 365)];
          }
          break;
        }
      }
    }


    LAND_leafAreaIndex_(_,g) = LAND_leafAreaIndexMin_(g) + LAND_leafAreaRatio_ * range_LAI(g);

  }

  return LAND_leafAreaIndex_;
}



//' @rdname process
//' @param param_INFILT_hbv_beta <0.001, 5> parameters for [infilt_HBV()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_HBV(
   NumericVector LAND_water_mm,
   NumericVector SOIL_water_mm,
   NumericVector SOIL_capacity_mm,
   NumericVector param_INFILT_hbv_beta
)
{
 NumericVector soil_diff_mm, infilt_water_mm, k_, limit_mm;

 soil_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
 limit_mm = ifelse(soil_diff_mm > LAND_water_mm, LAND_water_mm, soil_diff_mm);

 k_ = (1 - vecpow(SOIL_water_mm / SOIL_capacity_mm, param_INFILT_hbv_beta));

 infilt_water_mm = LAND_water_mm * k_;

 return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}


//' @rdname process
//' @param param_PERCOLA_arn_thresh <0.1, 0.9> coefficient parameter for [percola_ThreshPow()]
//' @param param_PERCOLA_arn_k <0.1, 1> exponential parameter for [percola_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector percola_Arno(
   NumericVector SOIL_water_mm,
   NumericVector SOIL_capacity_mm,
   NumericVector SOIL_potentialPercola_mm,
   NumericVector param_PERCOLA_arn_thresh,
   NumericVector param_PERCOLA_arn_k
)
{
 NumericVector percola_, percola_1, percola_2, Ws_Wc;
 Ws_Wc = SOIL_capacity_mm * param_PERCOLA_arn_thresh;
 percola_1 = param_PERCOLA_arn_k * SOIL_potentialPercola_mm / (SOIL_capacity_mm) * SOIL_water_mm;
 percola_2 = param_PERCOLA_arn_k * SOIL_potentialPercola_mm / (SOIL_capacity_mm) * SOIL_water_mm + SOIL_potentialPercola_mm * (1 - param_PERCOLA_arn_k) * pow((SOIL_water_mm - Ws_Wc) / (SOIL_capacity_mm - Ws_Wc),2);
 percola_ = ifelse(SOIL_water_mm < Ws_Wc, percola_1, percola_2);
 percola_ = ifelse(SOIL_potentialPercola_mm > Ws_Wc, SOIL_water_mm, percola_);
 percola_ = ifelse(percola_ > SOIL_potentialPercola_mm, SOIL_potentialPercola_mm, percola_);
 return ifelse(percola_ > SOIL_water_mm, SOIL_water_mm, percola_) ;
}

//' @rdname process
//' @param param_BASEFLOW_grf_gamma <2, 7> exponential parameter for [baseflow_GR4Jfix()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_GR4Jfix(
   NumericVector GROUND_water_mm,
   NumericVector GROUND_capacity_mm,
   NumericVector param_BASEFLOW_grf_gamma
)
{
 NumericVector baseflow_, k_;

 k_ = 1 - vecpow((1 + vecpow(GROUND_water_mm / GROUND_capacity_mm, param_BASEFLOW_grf_gamma)), -1.0 / param_BASEFLOW_grf_gamma);
 baseflow_ = k_ * GROUND_water_mm;

 return ifelse(baseflow_ > GROUND_water_mm, GROUND_water_mm, baseflow_) ;
}


//' @rdname process
//' @param param_BASEFLOW_sur_k <0.01, 1> coefficient parameter for [baseflow_SupplyRatio()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_SupplyRatio(
    NumericVector GROUND_water_mm,
    NumericVector param_BASEFLOW_sur_k
)
{

  return param_BASEFLOW_sur_k * GROUND_water_mm;

}


//' @rdname process
//' @param param_Lake_acp_storeFactor <uknow> parameter for [lake_AcceptPow()],
//' @param param_Lake_acp_gamma <uknow> parameter for [lake_AcceptPow()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector lake_AcceptPow(
   NumericVector Lake_water_m3,
   NumericVector Lake_inflow_m3,
   NumericVector Lake_capacity_m3,
   NumericVector param_Lake_acp_storeFactor,
   NumericVector param_Lake_acp_gamma
)
{
 Lake_water_m3 += Lake_inflow_m3;

 NumericVector Lake_outflow_m3 = (1 / param_Lake_acp_storeFactor) * vecpow(pmin(Lake_water_m3 / Lake_capacity_m3, 1), param_Lake_acp_gamma);
 Lake_water_m3 += -Lake_outflow_m3;

 NumericVector Lake_overflow_m3 = pmax(Lake_water_m3 -  Lake_capacity_m3, 0);
 Lake_water_m3 = pmin(Lake_water_m3, Lake_capacity_m3);

 Lake_outflow_m3 += Lake_overflow_m3;

 return (Lake_outflow_m3);
}


//' @rdname process
//' @export
// [[Rcpp::export]]
NumericVector river_LinearResorvoir(
   NumericVector RIVER_water_m3,
   NumericVector RIVER_inflow_m3,
   NumericVector RIVER_velocity_km,
   NumericVector RIVER_length_km
)
{

 NumericVector RIVER_paramK_TS = RIVER_length_km / RIVER_velocity_km; // pmax(RIVER_length_km / RIVER_velocity_km, 1.);

 // return RIVER_water_m3 * (1 / (RIVER_paramK_TS + 0.5)) + RIVER_inflow_m3 * (0.5 / (RIVER_paramK_TS + 0.5));
 return RIVER_water_m3 * (1 - exp(-1. / RIVER_paramK_TS)) + RIVER_inflow_m3 * (1 - RIVER_paramK_TS * (1 - exp(-1. / RIVER_paramK_TS)));
}

//' @rdname process
//' @param param_Riverlak_lin_storeFactor <uknow> parameter for [riverlak_LinearResorvoir()],
//' @export
// [[Rcpp::export]]
NumericVector riverlak_LinearResorvoir(
   NumericVector Riverlak_water_m3,
   NumericVector Riverlak_inflow_m3,
   NumericVector Riverlak_capacity_m3,
   NumericVector param_Riverlak_lin_storeFactor
)
{



 NumericVector Riverlak_outflow_m3 = Riverlak_water_m3 * (1 - exp(-1. / param_Riverlak_lin_storeFactor)) + Riverlak_inflow_m3 * (1 - param_Riverlak_lin_storeFactor * (1 - exp(-1. / param_Riverlak_lin_storeFactor)));


 NumericVector Riverlak_water_New = pmin(Riverlak_water_m3 + Riverlak_inflow_m3 - Riverlak_outflow_m3, Riverlak_capacity_m3);
 Riverlak_water_New = pmax(Riverlak_water_New, 0);
 Riverlak_outflow_m3 = Riverlak_water_m3 + Riverlak_inflow_m3 - Riverlak_water_New;

 return (Riverlak_outflow_m3);
}

//' @rdname process
//' @param param_Reservoi_han_alpha <uknow> parameter for [reservoi_Hanasaki()],
//' @param param_Reservoi_han_kDemand <uknow> parameter for [reservoi_Hanasaki()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector reservoi_Hanasaki(
   NumericVector Reservoi_water_m3,
   NumericVector Reservoi_inflow_m3,
   NumericVector Reservoi_capacity_m3,
   NumericVector Reservoi_demand_m3,
   NumericVector Reservoi_yearInflow_m3,
   NumericVector Reservoi_yearDemand_m3,
   NumericVector &Reservoi_yearRelase_m3,
   LogicalVector Reservoi_isOperateStart_01,
   LogicalVector Reservoi_isIrrigate_01,
   NumericVector param_Reservoi_han_alpha,
   NumericVector param_Reservoi_han_kDemand
)
{
 Reservoi_water_m3 += Reservoi_inflow_m3;

 NumericVector Reservoi_releaseCoefficient_ = Reservoi_water_m3 / (param_Reservoi_han_alpha * Reservoi_capacity_m3); // eq-3
 Reservoi_releaseCoefficient_ = pmax(Reservoi_releaseCoefficient_, 0.1);
 Reservoi_yearRelase_m3 = ifelse(Reservoi_isOperateStart_01, Reservoi_releaseCoefficient_ * Reservoi_yearInflow_m3, Reservoi_yearRelase_m3); // eq-2

 NumericVector Reservoi_releaseProvis_m3 = Reservoi_yearInflow_m3; // eq-4
 Reservoi_releaseProvis_m3 = ifelse(Reservoi_isIrrigate_01,
                                     0.5 * Reservoi_yearInflow_m3 * (1 + param_Reservoi_han_kDemand * Reservoi_demand_m3 / Reservoi_yearDemand_m3),
                                     Reservoi_releaseProvis_m3);  // eq-5
 Reservoi_releaseProvis_m3 = ifelse(Reservoi_yearDemand_m3 < 0.5 * Reservoi_yearInflow_m3,
                                     Reservoi_yearInflow_m3 + param_Reservoi_han_kDemand * Reservoi_demand_m3 - Reservoi_yearDemand_m3,
                                     Reservoi_releaseProvis_m3);  // eq-5
 NumericVector Reservoi_inflowRatio_ = Reservoi_capacity_m3 / Reservoi_yearInflow_m3; // eq-7
 NumericVector temp_inflowRatio_ = 4 * Reservoi_inflowRatio_ * Reservoi_inflowRatio_; // eq-7

 return ifelse(Reservoi_inflowRatio_ > 0.5,
               Reservoi_releaseCoefficient_ * Reservoi_releaseProvis_m3,
               temp_inflowRatio_ * Reservoi_releaseCoefficient_ * Reservoi_releaseProvis_m3 + (1 - temp_inflowRatio_) * Reservoi_inflow_m3); //  eq-7
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
   NumericVector CONFLUEN_cellInflow_m3,
   NumericVector &RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int
)
{

 int n_Cell = CONFLUEN_cellInflow_m3.size();
 NumericVector confluen_outflow_m3(n_Cell),
 step_RiverOutflow_m3, step_RiverlakeOutflow_m3;

 IntegerVector idx_Cell_Step,
 idx_RiverLake_Step, idx_Step_Riverlake;
 int n_Step = CELL_cellNumberStep_int.size();

 // Step i later with Inflow
 for (int i_Step = 0; i_Step < n_Step; i_Step++)
 {

   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   NumericVector step_UpstreamInflow_m3(idx_Cell_Step.size(), 0.);
   // Inflow upstream
   if (i_Step > 0) {

     step_UpstreamInflow_m3 = inflow_add(
       confluen_outflow_m3,
       CELL_inflowCellNumberStep_int[i_Step]
     );

   }


   // river segment
   NumericVector step_CellInflow = subset_get(CONFLUEN_cellInflow_m3, idx_Cell_Step),
     step_RiverWater = subset_get(RIVER_water_m3, idx_Cell_Step),
     step_RiverInflow = step_UpstreamInflow_m3 + step_CellInflow;

   step_RiverOutflow_m3 = river_LinearResorvoir(
     step_RiverWater,
     step_RiverInflow,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   NumericVector step_RIVER_Water_New = pmax(step_RiverWater + step_RiverInflow - step_RiverOutflow_m3, 0.0);
   NumericVector step_RIVER_Outflow_New = step_RiverWater + step_RiverInflow - step_RIVER_Water_New;
   subset_put(confluen_outflow_m3, idx_Cell_Step, step_RIVER_Outflow_New);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);


 }


 return confluen_outflow_m3;

}

//' @rdname process
//' @export
// [[Rcpp::export]]
NumericVector confluen_WaterGAP3_L(
   NumericVector CONFLUEN_cellInflow_m3,
   NumericVector &RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   IntegerVector Riverlak_cellNumber_int,
   NumericVector &Riverlak_water_m3,
   NumericVector Riverlak_capacity_m3,
   NumericVector param_Riverlak_lin_storeFactor
)
{

 int n_Cell = CONFLUEN_cellInflow_m3.size();
 NumericVector confluen_outflow_m3(n_Cell),
 step_RiverOutflow_m3, step_RiverlakeOutflow_m3;

 IntegerVector idx_Cell_Step,
 idx_RiverLake_Step, idx_Step_Riverlake;
 int n_Step = CELL_cellNumberStep_int.size();

 // Step i later with Inflow
 for (int i_Step = 0; i_Step < n_Step; i_Step++)
 {

   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   NumericVector step_UpstreamInflow_m3(idx_Cell_Step.size(), 0.);
   // Inflow upstream
   if (i_Step > 0) {

     step_UpstreamInflow_m3 = inflow_add(
       confluen_outflow_m3,
       CELL_inflowCellNumberStep_int[i_Step]
     );

   }


   // river segment
   NumericVector step_CellInflow = subset_get(CONFLUEN_cellInflow_m3, idx_Cell_Step),
     step_RiverWater = subset_get(RIVER_water_m3, idx_Cell_Step),
     step_RiverInflow = step_UpstreamInflow_m3 + step_CellInflow;

   step_RiverOutflow_m3 = river_LinearResorvoir(
     step_RiverWater,
     step_RiverInflow,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   NumericVector step_RIVER_Water_New = pmax(step_RiverWater + step_RiverInflow - step_RiverOutflow_m3, 0.0);
   NumericVector step_RIVER_Outflow_New = step_RiverWater + step_RiverInflow - step_RIVER_Water_New;
   subset_put(confluen_outflow_m3, idx_Cell_Step, step_RIVER_Outflow_New);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);

   // global lake (riverlake)
   idx_RiverLake_Step = get_idx_cell(Riverlak_cellNumber_int, idx_Cell_Step);
   idx_Step_Riverlake = get_idx_step(Riverlak_cellNumber_int, idx_Cell_Step);
   if (idx_RiverLake_Step.size() > 0) {
     NumericVector step_RiverlakeWater= subset_get(Riverlak_water_m3, idx_RiverLake_Step),
       step_RiverlakeInflow = subset_get(step_UpstreamInflow_m3 + step_CellInflow, idx_Step_Riverlake);


     step_RiverlakeOutflow_m3 = riverlak_LinearResorvoir(
       step_RiverlakeWater,
       step_RiverlakeInflow,
       subset_get(Riverlak_capacity_m3, idx_RiverLake_Step),
       subset_get(param_Riverlak_lin_storeFactor, idx_RiverLake_Step)
     );
     subset_put(confluen_outflow_m3, idx_RiverLake_Step, step_RiverlakeOutflow_m3);
     subset_put(Riverlak_water_m3, idx_RiverLake_Step,  step_RiverlakeWater + step_RiverlakeInflow - step_RiverlakeOutflow_m3);
   }



 }


 return confluen_outflow_m3;

}



//' @rdname process
//' @export
// [[Rcpp::export]]
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
)
{

 int n_Cell = CONFLUEN_cellInflow_m3.size();
 NumericVector confluen_outflow_m3(n_Cell), step_RiverOutflow_m3,
 step_RiverlakeOutflow_m3, step_ReservoirOutflow_m3;

 IntegerVector idx_Cell_Step,
 idx_Riverlake_Step, idx_Step_Riverlake,
 idx_Reservoir_Step, idx_Step_Reservoir;
 int n_Step = CELL_cellNumberStep_int.size();

 // Step i later with Inflow
 for (int i_Step = 0; i_Step < n_Step; i_Step++)
 {

   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   NumericVector step_UpstreamInflow_m3(idx_Cell_Step.size(), 0.);
   // Inflow upstream
   if (i_Step > 0) {

     step_UpstreamInflow_m3 = inflow_add(
       confluen_outflow_m3,
       CELL_inflowCellNumberStep_int[i_Step]
     );

   }


   // river segment
   NumericVector step_CellInflow = subset_get(CONFLUEN_cellInflow_m3, idx_Cell_Step),
     step_RiverWater= subset_get(RIVER_water_m3, idx_Cell_Step),
     step_RiverInflow = step_UpstreamInflow_m3 + step_CellInflow;

   step_RiverOutflow_m3 = river_LinearResorvoir(
     step_RiverWater,
     step_RiverInflow,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   NumericVector step_RiverInOut = pmax(step_RiverInflow - step_RiverOutflow_m3, 0.0);
   subset_put(confluen_outflow_m3, idx_Cell_Step, step_RiverOutflow_m3);
   subset_put(RIVER_water_m3, idx_Cell_Step, step_RiverWater + step_RiverInOut);

   // global lake (riverlake)
   idx_Riverlake_Step = get_idx_cell(Riverlak_cellNumber_int, idx_Cell_Step);
   idx_Step_Riverlake = get_idx_step(Riverlak_cellNumber_int, idx_Cell_Step);
   if (idx_Riverlake_Step.size() > 0) {
     NumericVector step_RiverlakeWater= subset_get(Riverlak_water_m3, idx_Riverlake_Step),
       step_RiverlakeInflow = subset_get(step_UpstreamInflow_m3 + step_CellInflow, idx_Step_Riverlake);


     step_RiverlakeOutflow_m3 = riverlak_LinearResorvoir(
       step_RiverlakeWater,
       step_RiverlakeInflow,
       subset_get(Riverlak_capacity_m3, idx_Riverlake_Step),
       subset_get(param_Riverlak_lin_storeFactor, idx_Riverlake_Step)
     );
     NumericVector step_RiverlakeInOut = pmax(step_RiverlakeInflow - step_RiverlakeOutflow_m3, 0.0);
     subset_put(confluen_outflow_m3, idx_Riverlake_Step, step_RiverlakeOutflow_m3);
     subset_put(Riverlak_water_m3, idx_Riverlake_Step, step_RiverlakeWater + step_RiverlakeInOut);
   }


   // Reservior
   idx_Reservoir_Step = get_idx_cell(Reservoi_cellNumber_int, idx_Cell_Step);
   idx_Step_Reservoir = get_idx_step(Reservoi_cellNumber_int, idx_Cell_Step);
   if (idx_Reservoir_Step.size() > 0 ) {
     NumericVector step_ReservoirYearRealse = subset_get(Reservoi_yearRelase_m3, idx_Reservoir_Step),
       step_ReservoirWater= subset_get(Reservoi_water_m3, idx_Reservoir_Step),
       step_ReservoirInflow = subset_get(step_UpstreamInflow_m3 + step_CellInflow, idx_Step_Reservoir);
     step_ReservoirOutflow_m3 = reservoi_Hanasaki(
       step_ReservoirWater,
       step_ReservoirInflow,
       subset_get(Reservoi_capacity_m3, idx_Reservoir_Step),
       subset_get(Reservoi_demand_m3, idx_Reservoir_Step),
       subset_get(Reservoi_yearInflow_m3, idx_Reservoir_Step),
       subset_get(Reservoi_yearDemand_m3, idx_Reservoir_Step),
       step_ReservoirYearRealse,
       subset_get_logical(Reservoi_isOperateStart_01, idx_Reservoir_Step),
       subset_get_logical(Reservoi_isIrrigate_01, idx_Reservoir_Step),
       subset_get(param_Reservoi_han_alpha, idx_Reservoir_Step),
       subset_get(param_Reservoi_han_kDemand, idx_Reservoir_Step)
     );
     NumericVector step_ReservoirInOut = pmax(step_ReservoirInflow - step_ReservoirOutflow_m3, 0.0);
     subset_put(Reservoi_yearRelase_m3, idx_Reservoir_Step, step_ReservoirYearRealse);
     subset_put(confluen_outflow_m3, idx_Reservoir_Step, step_ReservoirOutflow_m3);
     subset_put(Reservoi_water_m3, idx_Reservoir_Step, step_ReservoirWater + step_ReservoirInOut);
   }


 }


 return confluen_outflow_m3;

}






