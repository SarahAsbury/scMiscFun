
#' Apply Gating Conditions Based on Signs
#'
#' @description
#' Applies gating conditions to a SingleCellExperiment object based on specified
#' markers and gate signs. The function evaluates cells against boundary conditions
#' for two markers, determining which cells fall within the specified quadrant.
#'
#' @param sce A SingleCellExperiment object containing the data to be gated
#' @param assay_name Character string specifying which assay in the SingleCellExperiment to use
#' @param marker1 Character string specifying the name of the first marker (x-axis)
#' @param marker2 Character string specifying the name of the second marker (y-axis)
#' @param gate_signs Character string of length 2 specifying the gate directions:
#'        \itemize{
#'          \item First character (x-axis): '+' for greater than, '-' for less than
#'          \item Second character (y-axis): '+' for greater than, '-' for less than
#'        }
#' @param boundaries List containing gate boundary positions:
#'        \itemize{
#'          \item x: numeric value for x-axis boundary
#'          \item y: numeric value for y-axis boundary
#'        }
#'
#' @return A logical vector with TRUE for cells passing both gating conditions
#'
#' @details
#' The function evaluates two conditions based on the gate_signs parameter:
#' \itemize{
#'   \item For '+' signs, cells must be >= the boundary value
#'   \item For '-' signs, cells must be <= the boundary value
#' }
#' The final result is the logical AND of both conditions.
#'
#' @examples
#' \dontrun{
#' # Get cells in the upper-right quadrant (++):
#' passing_cells <- apply_gate_signs(
#'   sce = sce_object,
#'   assay_name = "protein",
#'   marker1 = "CD3",
#'   marker2 = "CD4",
#'   gate_signs = "++",
#'   boundaries = list(x = 2.5, y = 3.0)
#' )
#'
#' # Get cells in the upper-left quadrant (-+):
#' passing_cells <- apply_gate_signs(
#'   sce = sce_object,
#'   assay_name = "protein",
#'   marker1 = "CD3",
#'   marker2 = "CD4",
#'   gate_signs = "-+",
#'   boundaries = list(x = 2.5, y = 3.0)
#' )
#' }
#'
#' @importFrom SingleCellExperiment assay
#' @importFrom stringr str_sub
#'
#' @export



apply_gate_signs <- function(sce, assay_name, marker1, marker2, gate_signs, boundaries) {
  # X-axis condition
  x_condition <-
    if(str_sub(gate_signs, 1, 1) == "+") {
      assay(sce, assay_name)[marker1, ] >= boundaries$x
    } else if(str_sub(gate_signs, 1, 1) == "-") {
      assay(sce, assay_name)[marker1, ] <= boundaries$x
    } else {
      stop("Invalid x-axis gate sign")
    }

  # Y-axis condition
  y_condition <-
    if(str_sub(gate_signs, 2, 2) == "+") {
      assay(sce, assay_name)[marker2, ] >= boundaries$y
    } else if(str_sub(gate_signs, 2, 2) == "-") {
      assay(sce, assay_name)[marker2, ] <= boundaries$y
    } else {
      stop("Invalid y-axis gate sign")
    }

  # Combine conditions
  x_condition & y_condition
}
