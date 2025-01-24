// Defines a header file containing function for WATERGAP3_EvaCal/
#ifndef WATERGAP3_EvaCal
#define WATERGAP3_EvaCal

#include "00utilis.h"
double eva_NSE(NumericVector num_Sim, NumericVector num_Obs);
double eva_KGE(NumericVector num_Sim, NumericVector num_Obs,
               double factor_r = 1.0, double factor_alpha = 1.0, double factor_beta = 1.0);
List cali_DDS(Function fitness,
              List lst_OtherData,
              NumericVector x_Min, NumericVector x_Max,
              Nullable<NumericVector> x_Init = R_NilValue, int max_iter = 100,
              double r = 0.2);

#endif
