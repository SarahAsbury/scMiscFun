#' plot_doublets
#'
#' plot doublet scores. designed to work with scds doublet scores calculated in run_sc_qc pipeline.
#'
#' @param sce sce object
#' @param plot_by string. single variable to plot by
#' @param doublet_cols character vector. name of doublet score columns. default should work if using scdcs.
#' @return plot
#' @examples
#' plot_doublets(sce = sce, plot_by = params$sample_col)
#' @import dplyr
#' @import purrr
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @import scater
#' @import cowplot
#' @export

plot_doublets <- function(sce,
                          plot_by, # character string only
                          doublet_cols = c("cxds_score", "bcds_score", "hybrid_score")
){

  # individual scores
  score_plots <- map(doublet_cols,
                     ~sce %>%
                       plotColData(x = plot_by, y = .x, colour_by = plot_by) +
                       theme(axis.text.x = element_blank()) +
                       theme(legend.position = "none")
  ) %>%
    cowplot::plot_grid(plotlist = ., ncol = length(doublet_cols))

  # scoreA x scoreB
  combinations <- combn(doublet_cols, 2) %>% data.frame
  compare_scores <- map(combinations,

                        ~sce %>%
                          plotColData(x = .x[1], y = .x[2], colour_by = plot_by) +
                          theme(legend.position = "none")
  ) %>%
    cowplot::plot_grid(plotlist = ., ncol = ncol(combinations))

  # legend
  legend <- plotColData(sce, x = plot_by, y = doublet_cols[1], colour_by = plot_by) %>%
    cowplot::get_legend()

  # final plot
  out <- cowplot::plot_grid(score_plots,
                            compare_scores,
                            ncol = 1) %>%
    cowplot::plot_grid(legend, rel_widths = c(5,1))

  return(out)
}
