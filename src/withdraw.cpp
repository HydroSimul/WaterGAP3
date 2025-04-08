#include <Rcpp.h>
#include <algorithm> // for sort
using namespace Rcpp;

//' Withdraw water
//' @rdname withdraw
//'
//' @description
//' [withdraw_SingleCell] This function handles water withdrawal from a single cell's water storage.
//' It updates the withdrawal deficit and the remaining water volume.
//'
//' @param CELL_withdrawal_m3 Vector of withdrawal requirements in cubic meters.
//'   This represents water demand and will be updated to reflect remaining deficit.
//' @param CELL_water_m3 Vector of available water volumes in cubic meters.
//'   This will be updated to reflect remaining water after withdrawal.
//' @return None. Parameters are updated by reference.
//' @export
// [[Rcpp::export]]
void withdraw_SingleCell(NumericVector& CELL_withdrawal_m3,
                         NumericVector& CELL_water_m3) {
  NumericVector withdrawActual = ifelse(CELL_withdrawal_m3 > CELL_water_m3,
                                        CELL_water_m3,
                                        CELL_withdrawal_m3);

  // Update the deficit (passed by reference)
  CELL_withdrawal_m3 += -withdrawActual;

  // Updated cell water
  CELL_water_m3 += -withdrawActual;
}

//' [withdrawSurface_AroundMax] This function identifies the surrounding cell with the maximum water availability
//' and withdraws water from it to satisfy the demand.
//' @name withdraw
//' @param CELL_withdrawal_m3 Vector of water withdrawal requirements in cubic meters.
//'   Will be updated to reflect remaining deficits after withdrawal.
//' @param RIVER_water_m3 Vector of available river water volumes in cubic meters.
//'   Will be updated to reflect remaining water after withdrawal.
//' @param RESERVOIR_water_m3 Vector of available reservoir water volumes in cubic meters.
//'   Will be updated to reflect remaining water after withdrawal.
//' @param RIVERLAK_water_m3 Vector of available river-lake water volumes in cubic meters.
//'   Will be updated to reflect remaining water after withdrawal.
//' @param LAKE_water_m3 Vector of available lake water volumes in cubic meters.
//'   Will be updated to reflect remaining water after withdrawal.
//' @param CELL_cellNumberAround_int Matrix of cell indices that define surrounding cells.
//'   Each column corresponds to a cell, with rows containing indices of surrounding cells.
//' @return None. Parameters are updated by reference.
//' @export
// [[Rcpp::export]]
void withdrawSurface_AroundMax(NumericVector& CELL_withdrawal_m3,
                               NumericVector& RIVER_water_m3,
                               NumericVector& RESERVOIR_water_m3,
                               NumericVector& RIVERLAK_water_m3,
                               NumericVector& LAKE_water_m3,
                               IntegerMatrix CELL_cellNumberAround_int) {

  int n_Spat = CELL_withdrawal_m3.size();

  for(int i_Spat = 0; i_Spat < n_Spat; i_Spat++) {

    if(CELL_withdrawal_m3[i_Spat] <= 0) continue;

    // Get surrounding cell indices
    IntegerVector idx_Around = CELL_cellNumberAround_int(_, i_Spat) - 1; // Convert to 0-based index
    idx_Around = idx_Around[idx_Around >= 0];
    int n_Around = idx_Around.size();
    if(n_Around == 0) continue;

    // Calculate total water in surrounding cells
    NumericVector allWater_Around(n_Around);
    for(int j = 0; j < n_Around; j++) {
      int cell_idx = idx_Around[j];
      allWater_Around[j] = RIVER_water_m3[cell_idx] +
        RESERVOIR_water_m3[cell_idx] +
        RIVERLAK_water_m3[cell_idx] +
        LAKE_water_m3[cell_idx];
    }

    // Find cell with maximum water
    int idx_MaxAround = 0;
    double max_water = allWater_Around[0];
    for(int j = 1; j < n_Around; j++) {
      if(allWater_Around[j] > max_water) {
        max_water = allWater_Around[j];
        idx_MaxAround = j;
      }
    }

    int target_cell = idx_Around[idx_MaxAround]; // 0-based index
    double withdrawal = CELL_withdrawal_m3[i_Spat];
    double river = RIVER_water_m3[target_cell];
    double reservoir = RESERVOIR_water_m3[target_cell];
    double riverlak = RIVERLAK_water_m3[target_cell];
    double lake = LAKE_water_m3[target_cell];
    double total_water = river + reservoir + riverlak + lake;

    // Apply withdrawal rules
    if(withdrawal > total_water) {
      RIVER_water_m3[target_cell] = 0;
      RESERVOIR_water_m3[target_cell] = 0;
      RIVERLAK_water_m3[target_cell] = 0;
      LAKE_water_m3[target_cell] = 0;
      CELL_withdrawal_m3[i_Spat] -= total_water;
    }
    else if(withdrawal > (total_water - lake)) {
      RIVER_water_m3[target_cell] = 0;
      RESERVOIR_water_m3[target_cell] = 0;
      RIVERLAK_water_m3[target_cell] = 0;
      LAKE_water_m3[target_cell] += (river + reservoir + riverlak - withdrawal);
      CELL_withdrawal_m3[i_Spat] = 0;
    }
    else if(withdrawal > (total_water - lake - riverlak)) {
      RIVER_water_m3[target_cell] = 0;
      RESERVOIR_water_m3[target_cell] = 0;
      RIVERLAK_water_m3[target_cell] += (river + reservoir - withdrawal);
      CELL_withdrawal_m3[i_Spat] = 0;
    }
    else if(withdrawal > (river + reservoir)) {
      RIVER_water_m3[target_cell] = 0;
      RESERVOIR_water_m3[target_cell] += (river - withdrawal);
      CELL_withdrawal_m3[i_Spat] = 0;
    }
    else {
      RIVER_water_m3[target_cell] = river - withdrawal;
      CELL_withdrawal_m3[i_Spat] = 0;
    }
  }
}

//' [withdrawSurface_Around] This function withdraws water from all surrounding cells proportionally
//' based on their water availability.
//' @name withdraw
//' @export
// [[Rcpp::export]]
void withdrawSurface_Around(NumericVector& CELL_withdrawal_m3,
                            NumericVector& RIVER_water_m3,
                            NumericVector& RESERVOIR_water_m3,
                            NumericVector& RIVERLAK_water_m3,
                            NumericVector& LAKE_water_m3,
                            IntegerMatrix CELL_cellNumberAround_int) {

  int n_Spat = CELL_withdrawal_m3.size();

  for(int i_Spat = 0; i_Spat < n_Spat; i_Spat++) {
    // Skip if there's no withdrawal needed
    if(CELL_withdrawal_m3[i_Spat] <= 0) continue;

    // Get surrounding cell indices (convert to 0-based)
    IntegerVector idx_Around = CELL_cellNumberAround_int(_, i_Spat) - 1;
    idx_Around = idx_Around[idx_Around >= 0];
    int n_Around = idx_Around.size();
    if (n_Around == 0) continue;

    // Skip cells with invalid indices
    IntegerVector valid_indices;
    for(int j = 0; j < n_Around; j++) {
      if(idx_Around[j] >= 0) {
        valid_indices.push_back(idx_Around[j]);
      }
    }

    // Skip if no valid cells
    if(valid_indices.size() == 0) continue;

    // Calculate total water in surrounding cells
    double total_available = 0;
    NumericVector totalWater_Around(valid_indices.size());

    for(int j = 0; j < valid_indices.size(); j++) {
      int idx = valid_indices[j];
      double water = RIVER_water_m3[idx] + RESERVOIR_water_m3[idx] +
        RIVERLAK_water_m3[idx] + LAKE_water_m3[idx];
      totalWater_Around[j] = water;
      total_available += water;
    }

    double withdrawal = CELL_withdrawal_m3[i_Spat];

    if(withdrawal > total_available) {
      // Case 1: Withdrawal exceeds total available water
      for(int j = 0; j < valid_indices.size(); j++) {
        int idx = valid_indices[j];
        RIVER_water_m3[idx] = 0;
        RESERVOIR_water_m3[idx] = 0;
        RIVERLAK_water_m3[idx] = 0;
        LAKE_water_m3[idx] = 0;
      }
      CELL_withdrawal_m3[i_Spat] -= total_available;
    } else {
      // Case 2: Withdrawal can be satisfied from available water
      NumericVector weights = totalWater_Around / total_available;

      for(int j = 0; j < valid_indices.size(); j++) {
        int idx = valid_indices[j];
        double cell_withdrawal = withdrawal * weights[j];

        // Apply withdrawal prioritizing different water sources
        double river = RIVER_water_m3[idx];
        double reservoir = RESERVOIR_water_m3[idx];
        double riverlak = RIVERLAK_water_m3[idx];
        double lake = LAKE_water_m3[idx];
        double total_cell_water = river + reservoir + riverlak + lake;

        if(cell_withdrawal >= total_cell_water) {
          RIVER_water_m3[idx] = 0;
          RESERVOIR_water_m3[idx] = 0;
          RIVERLAK_water_m3[idx] = 0;
          LAKE_water_m3[idx] = 0;
        } else if(cell_withdrawal > river) {
          RIVER_water_m3[idx] = 0;
          double remaining = cell_withdrawal - river;

          if(remaining > reservoir) {
            RESERVOIR_water_m3[idx] = 0;
            remaining -= reservoir;

            if(remaining > riverlak) {
              RIVERLAK_water_m3[idx] = 0;
              LAKE_water_m3[idx] -= (remaining - riverlak);
            } else {
              RIVERLAK_water_m3[idx] -= remaining;
            }
          } else {
            RESERVOIR_water_m3[idx] -= remaining;
          }
        } else {
          RIVER_water_m3[idx] -= cell_withdrawal;
        }
      }

      CELL_withdrawal_m3[i_Spat] = 0;
    }
  }
}

//' [withdrawSurface_WithdrawNet] This function withdraws water by following a predefined network of cells
//' from which water can be withdrawn.
//' @name withdraw
//' @param CELL_cellNumberWithdrawNet_int Matrix defining the withdrawal network.
//'   Each column corresponds to a cell, with rows containing indices of cells
//'   from which water can be withdrawn in order of priority.
//' @export
// [[Rcpp::export]]
void withdrawSurface_WithdrawNet(NumericVector& CELL_withdrawal_m3,
                                 NumericVector& RIVER_water_m3,
                                 NumericVector& RESERVOIR_water_m3,
                                 NumericVector& RIVERLAK_water_m3,
                                 NumericVector& LAKE_water_m3,
                                 IntegerMatrix CELL_cellNumberWithdrawNet_int) {

  int n_Spat = CELL_withdrawal_m3.size();

  for(int i_Spat = 0; i_Spat < n_Spat; i_Spat++) {
    // Skip if there's no withdrawal needed
    if(CELL_withdrawal_m3[i_Spat] <= 0) continue;

    // Get withdrawal network indices (convert to 0-based)
    IntegerVector idx_WithdrawNet = CELL_cellNumberWithdrawNet_int(_, i_Spat) - 1;

    int idx_WD = 0;
    while(CELL_withdrawal_m3[i_Spat] > 0 && idx_WD < idx_WithdrawNet.size()) {
      int idx_Spat_WD = idx_WithdrawNet[idx_WD];

      // Skip if index is invalid (0 in R becomes -1 in C++ after conversion)
      if(idx_Spat_WD < 0) {
        idx_WD++;
        continue;
      }

      double river = RIVER_water_m3[idx_Spat_WD];
      double reservoir = RESERVOIR_water_m3[idx_Spat_WD];
      double riverlak = RIVERLAK_water_m3[idx_Spat_WD];
      double lake = LAKE_water_m3[idx_Spat_WD];
      double total_Water = river + reservoir + riverlak + lake;
      double withdrawal = CELL_withdrawal_m3[i_Spat];

      if(withdrawal > total_Water) {
        // Case 1: Withdrawal exceeds total available water
        RIVER_water_m3[idx_Spat_WD] = 0;
        RESERVOIR_water_m3[idx_Spat_WD] = 0;
        RIVERLAK_water_m3[idx_Spat_WD] = 0;
        LAKE_water_m3[idx_Spat_WD] = 0;
        CELL_withdrawal_m3[i_Spat] -= total_Water;
      }
      else if(withdrawal > (total_Water - lake)) {
        // Case 2: Withdrawal affects lake storage
        RIVER_water_m3[idx_Spat_WD] = 0;
        RESERVOIR_water_m3[idx_Spat_WD] = 0;
        RIVERLAK_water_m3[idx_Spat_WD] = 0;
        LAKE_water_m3[idx_Spat_WD] = lake - (withdrawal - (river + reservoir + riverlak));
        CELL_withdrawal_m3[i_Spat] = 0;
        break;
      }
      else if(withdrawal > (total_Water - lake - riverlak)) {
        // Case 3: Withdrawal affects river-lake storage
        RIVER_water_m3[idx_Spat_WD] = 0;
        RESERVOIR_water_m3[idx_Spat_WD] = 0;
        RIVERLAK_water_m3[idx_Spat_WD] = riverlak - (withdrawal - (river + reservoir));
        CELL_withdrawal_m3[i_Spat] = 0;
        break;
      }
      else if(withdrawal > river) {
        // Case 4: Withdrawal affects reservoir storage
        RIVER_water_m3[idx_Spat_WD] = 0;
        RESERVOIR_water_m3[idx_Spat_WD] = reservoir - (withdrawal - river);
        CELL_withdrawal_m3[i_Spat] = 0;
        break;
      }
      else {
        // Case 5: Withdrawal can be satisfied from river only
        RIVER_water_m3[idx_Spat_WD] = river - withdrawal;
        CELL_withdrawal_m3[i_Spat] = 0;
        break;
      }

      idx_WD++;
    }
  }
}
