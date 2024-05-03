#' plot_manual_pairwise
#'
#' pairwise plots with manual gates. include reports of percentage of cells in each gated population.
#'
#' @param sce sce object
#' @param features characteor vector of the 2 features to plot
#' @param thresholds character vectors of the thresholds for each of the 2 features in the same order as features
#' @param exprs_values sce object assay to use
#' @param altExp_name sce object alt exp name to use
#' @return list of sce objects
#' @examples
#' threshold$nk$lognorm <- c(CD3 = 0.25, CD56 = 0)
#' gating$nk$lognorm <- plot_manual_pairwise(sce,
#'                                           features = c("CD3", "CD56"),
#'                                           thresholds = threshold$nk$lognorm,
#'                                           exprs_value = "logcounts")
#' @import dplyr
#' @import purrr
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @import scater
#' @export


plot_manual_pairwise <- function(sce,
                                 features,
                                 thresholds,
                                 altExp_name = "PROTEIN",
                                 exprs_value = "logcounts"
){

  # Extract feature names from features vector
  pos1 <- paste0(features[[1]], "_pos")
  pos2 <- paste0(features[[2]], "_pos")

  # Assign population by threshold
  cell_pop <- sce %>%
    altExp(altExp_name) %>%
    assay(exprs_value) %>%
    .[c(features[[1]], features[[2]]), ] %>%
    as.matrix() %>%
    t() %>%
    data.frame(check.names = F) %>%
    mutate(!!pos1 := ifelse(.data[[features[[1]]]] > thresholds[[1]], TRUE, FALSE),
           !!pos2 := ifelse(.data[[features[[2]]]] > thresholds[[2]], TRUE, FALSE)) %>%
    mutate(population = case_when(
      !!sym(pos1) & !!sym(pos2) ~ paste0(features[[1]], "+", features[[2]], "+"),
      !!sym(pos1) & !!sym(pos2) == F ~ paste0(features[[1]], "+"),
      !!sym(pos1) ==F & !!sym(pos2)  ~ paste0(features[[2]], "+"),
      TRUE ~ paste0(features[[1]], "-", features[[2]], "-")
    ))

  # population statistics (count + freq)
  pop_tab <- cell_pop %>% group_by(population) %>%
    tally %>%
    mutate(freq = n/sum(n))

  # format and print as percentage
  pop_tab %>%
    mutate(percent = paste0(round(freq * 100, digit = 2), "%")) %>%
    select(-n, -freq) %>%
    print()

  # plot
  # deciles <-  pop_tab %>% select(CD4, CD8) %>% map(~.x %>% quantile(probs =  seq(0, 1, by = 0.1)))
  sce %>% visualiseExprs(altExp_name = "PROTEIN",
                         exprs_value = exprs_value,
                         plot = "pairwise",
                         feature_subset = c(features),
                         threshold = thresholds)

  return(list(all = cell_pop, summary = pop_tab))

}



