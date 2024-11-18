#include "WaterGAP3.h"

//' WaterGAP3
//' @name WaterGAP3_H
//' @inheritParams all_vari
//' @inheritParams process
//' @return streamflow m3
//' @export
// [[Rcpp::export]]
List WaterGAP3_H(
   int n_time,
   int n_spat,
   NumericMatrix AtmoS_precipitation_mm,
   NumericMatrix AtmoS_temperature_Cel,
   NumericMatrix AtmoS_solarRadiat_MJ,
   NumericVector snoW_ice_mm,
   NumericVector soiL_water_mm,
   NumericVector soiL_capacity_mm,
   NumericVector soiL_potentialPercola_mm,
   NumericVector grounD_water_mm,
   NumericVector grounD_capacity_mm,
   NumericVector riveR_water_m3,
   NumericVector riveR_length_km,
   NumericVector riveR_velocity_km,
   NumericVector celL_landArea_km2,
   List celL_cellNumberStep_int,
   List celL_inflowCellNumberStep_int,
   NumericVector param_atmos_thr_Ts,
   NumericVector param_snow_fac_f,
   NumericVector param_snow_fac_Tmelt,
   NumericVector param_evatrans_tur_k,
   NumericVector param_evatrans_ubc_gamma,
   NumericVector param_infilt_hbv_beta,
   NumericVector param_percola_arn_k,
   NumericVector param_percola_arn_thresh,
   NumericVector param_baseflow_grf_gamma,
   bool if_allVariExport = false
)
{

 NumericVector AtmoS_potentialEvatrans_mm, AtmoS_rain_mm, AtmoS_snoW_mm,
               snoW_melt_mm, lanD_water_mm,
               soiL_evatranS_mm, soiL_infilt_mm, soiL_percola_mm,
               celL_outflow_m3;
 NumericMatrix LanD_runoff_mm(n_time, n_spat), GrounD_basefloW_mm(n_time, n_spat), RiveR_outflow_m3(n_time, n_spat);

 NumericMatrix OuT_snow(n_time, n_spat), OuT_evatransPot(n_time, n_spat), OuT_evatrans(n_time, n_spat),
 OuT_soilwater(n_time, n_spat), OuT_groundwater(n_time, n_spat),
 OuT_snowice(n_time, n_spat), OuT_snowmelt(n_time, n_spat),
 OuT_cellOutflow(n_time, n_spat),
 OuT_riverwater(n_time, n_spat);



 for (int i = 0; i < n_time; i++) {

   // Land
   // // snow
   AtmoS_snoW_mm = atmosSnow_ThresholdT(AtmoS_precipitation_mm(i, _), AtmoS_temperature_Cel(i, _), param_atmos_thr_Ts);
   AtmoS_rain_mm = AtmoS_precipitation_mm(i, _) - AtmoS_snoW_mm;

   // // PET
   AtmoS_potentialEvatrans_mm = evatransPotential_TurcWendling(AtmoS_temperature_Cel(i, _), AtmoS_solarRadiat_MJ(i, _), param_evatrans_tur_k);

   // // soil
   soiL_evatranS_mm = evatransActual_UBC(AtmoS_potentialEvatrans_mm, soiL_water_mm, soiL_capacity_mm, param_evatrans_ubc_gamma);
   soiL_water_mm += - soiL_evatranS_mm;
   lanD_water_mm = AtmoS_rain_mm;

   // // Snow melt
   snoW_melt_mm = snowMelt_Factor(snoW_ice_mm, AtmoS_temperature_Cel(i, _), param_snow_fac_f, param_snow_fac_Tmelt);
   lanD_water_mm += snoW_melt_mm;
   snoW_ice_mm += -snoW_melt_mm;
   snoW_ice_mm += AtmoS_snoW_mm;

   // // soil infiltration
   soiL_infilt_mm = infilt_HBV(lanD_water_mm, soiL_water_mm, soiL_capacity_mm, param_infilt_hbv_beta);
   soiL_water_mm += soiL_infilt_mm;
   LanD_runoff_mm(i, _) = lanD_water_mm - soiL_infilt_mm;

   // // soil percolation
   soiL_percola_mm = percola_Arno(soiL_water_mm, soiL_capacity_mm, soiL_potentialPercola_mm, param_percola_arn_thresh, param_percola_arn_k);
   grounD_water_mm += soiL_percola_mm;
   soiL_water_mm += - soiL_percola_mm;

   // // baseflow
   NumericVector basefloW_temp = ifelse(grounD_water_mm < grounD_capacity_mm, 0, grounD_water_mm - grounD_capacity_mm);

   // // ground water
   grounD_water_mm = ifelse(grounD_water_mm < grounD_capacity_mm,grounD_water_mm, grounD_capacity_mm);
   GrounD_basefloW_mm(i, _) = baseflow_GR4Jfix(grounD_water_mm, grounD_capacity_mm, param_baseflow_grf_gamma);
   grounD_water_mm += - GrounD_basefloW_mm(i, _);
   GrounD_basefloW_mm(i, _) = GrounD_basefloW_mm(i, _) + basefloW_temp;


   celL_outflow_m3 = (LanD_runoff_mm(i, _) + GrounD_basefloW_mm(i, _)) * celL_landArea_km2 * 1000;



   RiveR_outflow_m3(i, _) = confluen_WaterGAP3(
     celL_outflow_m3,
     riveR_water_m3,
     riveR_length_km,
     riveR_velocity_km,
     celL_cellNumberStep_int,
     celL_inflowCellNumberStep_int
   );

   if (if_allVariExport) {
     OuT_cellOutflow(i, _) = celL_outflow_m3;
     OuT_snow(i, _) = AtmoS_snoW_mm;
     OuT_evatransPot(i, _) = AtmoS_potentialEvatrans_mm;
     OuT_evatrans(i, _) = soiL_evatranS_mm;
     OuT_soilwater(i, _) = soiL_water_mm;
     OuT_groundwater(i, _) = grounD_water_mm;
     OuT_snowice(i, _) = snoW_ice_mm;
     OuT_snowmelt(i, _) = snoW_melt_mm;
     OuT_riverwater(i, _) = riveR_water_m3;
   }

 }

 if (if_allVariExport) {
   return List::create(
     _["snowFall_mm"] = OuT_snow,
     _["evatransPot_mm"] = OuT_evatransPot,
     _["evatranS_mm"] = OuT_evatrans,
     _["soilwater_mm"] = OuT_soilwater,
     _["groundwater_mm"] = OuT_groundwater,
     _["snowice_mm"] = OuT_snowice,
     _["snowmelt_mm"] = OuT_snowmelt,
     _["runoff_mm"] = LanD_runoff_mm,
     _["basefloW_mm"] = GrounD_basefloW_mm,
     _["runoff_m3"] = OuT_cellOutflow,
     _["river_water_m3"] = OuT_riverwater,
     _["streamflow_m3"] = RiveR_outflow_m3
   );
 } else{
   return List::create(
     _["streamflow_m3"] = RiveR_outflow_m3
   );
 }


}
