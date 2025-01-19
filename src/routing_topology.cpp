#include <Rcpp.h>
#include <unordered_map>
#include <algorithm>
using namespace Rcpp;


//' @rdname routingtopology
//' @name get_inflow_cells
//' @title Get Inflow Cells
//' @description This function calculates inflow cells based on the outflow vector.
//' @param int_Outflow (vector of int) The cell number of the next cell. The cell number must range from 1 to the length of the cells.
//' If the cell has no outflow, the number should be set to itself.
//' @return A list where each element is an IntegerVector containing the inflow cells for each cell.
//' @export
// [[Rcpp::export]]
List get_inflow_cells(IntegerVector int_Outflow) {
  int n_Cell = int_Outflow.size();
  List lst_Inflow(n_Cell);

  // Initialize each element of lst_Inflow with an empty IntegerVector
  for (int i = 0; i < n_Cell; i++) {
    lst_Inflow[i] = IntegerVector();
  }

  for (int i = 0; i < n_Cell; i++) {
    IntegerVector inflow_i = as<IntegerVector>(lst_Inflow[i]);
    inflow_i.push_back(i + 1); // Start with the current cell (1-based index)
    lst_Inflow[i] = inflow_i;

    int num_Ori = i + 1;  // Start with the current cell (1-based index)
    int num_Next = int_Outflow[i]; // Get the next cell from the flow direction vector

    while (num_Next != num_Ori) {
      // Append the current cell to the inflow list of num_Next
      IntegerVector inflow_next = as<IntegerVector>(lst_Inflow[num_Next - 1]);
      inflow_next.push_back(i + 1);  // Adjust for 1-based indexing
      lst_Inflow[num_Next - 1] = inflow_next;

      // Update num_Ori and num_Next for the next iteration
      num_Ori = num_Next;
      num_Next = int_Outflow[num_Next - 1]; // Move to the next cell
    }
  }

  return lst_Inflow;
}

//' @rdname routingtopology
//' @name get_inflow_lastcell
//' @title Get Inflow Last Cell Matrix
//' @description This function creates a matrix of inflow cells for each cell based on the outflow vector.
//' @return A NumericMatrix where each row corresponds to a cell, and each column represents an inflow cell.
//' @export
// [[Rcpp::export]]
NumericMatrix get_inflow_lastcell(IntegerVector int_Outflow) {
  int n_Cell = int_Outflow.size();

  // Create a list to store inflow cells
  std::vector<std::vector<int>> lst_Inflow_LastCell(n_Cell);

  // Populate lst_Inflow_LastCell with indices
  for (int i = 0; i < n_Cell; i++) {
    for (int j = 0; j < n_Cell; j++) {
      if (int_Outflow[j] == i + 1) {  // Adjust for 1-based indexing
        lst_Inflow_LastCell[i].push_back(j + 1);  // Store 1-based index
      }
    }
  }

  // Find the maximum number of inflow cells for any cell
  int max_n_LastCell = 0;
  for (int i = 0; i < n_Cell; i++) {
    if (lst_Inflow_LastCell[i].size() > max_n_LastCell) {
      max_n_LastCell = lst_Inflow_LastCell[i].size();
    }
  }

  // Initialize the matrix with NA values
  NumericMatrix mat_Inflow_LastCell(n_Cell, max_n_LastCell);
  std::fill(mat_Inflow_LastCell.begin(), mat_Inflow_LastCell.end(), NA_REAL);

  // Fill the matrix with the inflow cells
  for (int i = 0; i < n_Cell; i++) {
    for (int j = 0; j < lst_Inflow_LastCell[i].size(); j++) {
      mat_Inflow_LastCell(i, j) = lst_Inflow_LastCell[i][j];
    }
  }

  return mat_Inflow_LastCell;
}


List get_step_cells(List lst_Flow_Cells) {
  int n_Cell = lst_Flow_Cells.size();

  // Step 1: Calculate the lengths of each element in lst_Flow_Cells
  IntegerVector length_Inflow(n_Cell);
  for (int i = 0; i < n_Cell; i++) {
    length_Inflow[i] = Rf_length(lst_Flow_Cells[i]);
  }

  // Step 2: Create a table to count occurrences of each length
  std::unordered_map<int, int> table_Count;
  for (int i = 0; i < n_Cell; i++) {
    table_Count[length_Inflow[i]]++;
  }

  // Step 3: Sort the lengths and assign step numbers
  std::vector<int> sorted_lengths(n_Cell);
  std::copy(length_Inflow.begin(), length_Inflow.end(), sorted_lengths.begin());
  std::sort(sorted_lengths.begin(), sorted_lengths.end());

  IntegerVector i_Step(n_Cell);
  int step = 1;
  for (int i = 0; i < n_Cell; i++) {
    if (i == 0 || sorted_lengths[i] != sorted_lengths[i - 1]) {
      i_Step[i] = step++;
    } else {
      i_Step[i] = i_Step[i - 1];
    }
  }

  // Step 4: Match length_Inflow to sorted_lengths and get corresponding step numbers
  IntegerVector idx_Step_Cell(n_Cell);
  for (int i = 0; i < n_Cell; i++) {
    auto it = std::find(sorted_lengths.begin(), sorted_lengths.end(), length_Inflow[i]);
    idx_Step_Cell[i] = i_Step[std::distance(sorted_lengths.begin(), it)];
  }

  // Step 5: Group cells by their step numbers
  List lst_Step_Cell(step - 1);
  for (int i = 0; i < n_Cell; i++) {
    int idx = idx_Step_Cell[i] - 1;  // Adjust for 0-based indexing
    if (lst_Step_Cell[idx] == R_NilValue) {
      lst_Step_Cell[idx] = IntegerVector::create(i + 1);  // Initialize with 1-based index
    } else {
      IntegerVector temp = as<IntegerVector>(lst_Step_Cell[idx]);
      temp.push_back(i + 1);  // Add 1-based index
      lst_Step_Cell[idx] = temp;
    }
  }

  return lst_Step_Cell;
}


List get_step_lastcell(List lst_Step_Cell, NumericMatrix mat_Inflow_LastCell) {
  int n_Step = lst_Step_Cell.size();
  List lst_Step_InflowLastCell(n_Step);

  // Set the first element to NA (equivalent in Rcpp is R_NilValue)
  lst_Step_InflowLastCell[0] = R_NilValue;

  for (int i = 1; i < n_Step; i++) {
    IntegerVector indices = lst_Step_Cell[i];
    int n_Indices = indices.size();
    NumericMatrix submatrix(n_Indices, mat_Inflow_LastCell.ncol());

    for (int j = 0; j < n_Indices; j++) {
      int idx = indices[j] - 1;  // Adjust for 0-based indexing in C++
      submatrix(j, _) = mat_Inflow_LastCell(idx, _);
    }

    lst_Step_InflowLastCell[i] = submatrix;
  }

  return lst_Step_InflowLastCell;
}


//' @rdname routingtopology
//' @name get_step_param
//' @title Get Step Parameters
//' @description This function returns a list of step cells and the corresponding last cell matrices.
//' @return A list containing the step cells and last cell matrices.
//' @export
// [[Rcpp::export]]
List get_step_param(IntegerVector int_Outflow) {
  // Step 1: Get inflow cells
  List lst_Flow_Cells = get_inflow_cells(int_Outflow);

  // Step 2: Get step cells based on inflow cells
  List lst_Step_Cell = get_step_cells(lst_Flow_Cells);

  // Step 3: Get inflow last cell matrix based on outflow data
  NumericMatrix mat_Inflow_LastCell = get_inflow_lastcell(int_Outflow);

  // Step 4: Get step last cell list based on step cells and inflow last cell matrix
  List lst_Step_LastCell = get_step_lastcell(lst_Step_Cell, mat_Inflow_LastCell);

  // Return the results as a list
  return List::create(
    Named("int_Cell") = lst_Step_Cell,
    Named("mat_LastCell") = lst_Step_LastCell
  );
}
