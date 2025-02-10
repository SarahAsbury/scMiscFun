#' Create Gate Visualization
#'
#' @param sce SingleCellExperiment object
#' @param marker1 First marker name (x-axis)
#' @param marker2 Second marker name (y-axis)
#' @param x_boundary Numeric value for vertical gate position
#' @param y_boundary Numeric value for horizontal gate position
#' @param population_name Character string for plot title
#' @param assay_name Character string specifying which assay to use
#'
#' @import ggplot2
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @import magrittr
#' @import dplyr
#' @import purrr

#'
#' @return ggplot object
#' @export
create_gate_plot <- function(sce,
                             marker1, marker2,
                             x_boundary, y_boundary,
                             population_name, assay_name
                             ) {
  plot_data <- data.frame(
    x = assay(sce, assay_name)[marker1, ],
    y = assay(sce, assay_name)[marker2, ]
  )

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.2) +
    geom_vline(xintercept = x_boundary, color = "red", linetype = "dashed") +
    geom_hline(yintercept = y_boundary, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = population_name,
      x = marker1,
      y = marker2
    )

  return(p)
}
