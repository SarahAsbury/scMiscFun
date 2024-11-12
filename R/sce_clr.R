#' sce_clr
#' Apply centered log-ratio transformation to SingleCellExperiment object
#'
#' Performs centered log-ratio (CLR) transformation on count data stored in an
#' alternative experiment of a SingleCellExperiment object. The function adds
#' pseudocounts before transformation and performs quality checks on the results.
#'
#' @param x A SingleCellExperiment object containing the data to transform
#' @param altExp_name Character string specifying the name of the alternative
#'   experiment containing the count data (default: "PROTEIN")
#'
#' @return A modified SingleCellExperiment object with:
#'   \itemize{
#'     \item Original counts stored in "rawcounts" assay
#'     \item CLR-transformed data stored in both "counts" and "logcounts" assays
#'   }
#'
#' @details
#' The function performs the following steps:
#' 1. Adds pseudocounts using the pseudo_count function
#' 2. Applies CLR transformation using the compositions package
#' 3. Performs quality checks on row and column names
#' 4. Stores the results in appropriate assay slots
#'
#' @import SingleCellExperiment
#' @importFrom compositions clr
#' @importFrom purrr walk2
#' @importFrom stringr str_replace
#' @importFrom magrittr %>%
#'
#' @examples
#' # Assuming 'sce' is a SingleCellExperiment object with protein data
#' sce <- sce_clr(sce, altExp_name = "PROTEIN")
#'
#' @export
sce_clr <- function(x, altExp_name = "PROTEIN") {
  psuedo_clr <- x %>%
    altExp(altExp_name) %>%
    counts %>%
    apply(MARGIN = 2, pseudo_count) #add psuedo counts
  calc_clr <- psuedo_clr %>% apply(MARGIN = 2, compositions::clr)

  # sanity check
  qc <- list(
    rows_equal = all(rownames(calc_clr) == rownames(x %>% altExp(altExp_name))),
    cols_equal = all((colnames(calc_clr) %>% str_replace("\\.", "-")) == colnames(x %>% altExp(altExp_name)))
  )

  walk2(qc, names(qc),
        ~if(!.x) stop(sprintf("%s qc check failed", .y))
  )

  # save raw
  assay(altExp(x, altExp_name), "rawcounts") <- x %>% altExp %>% counts

  # add clr to both assay slots
  assay(altExp(x, altExp_name), "counts") <- calc_clr
  assay(altExp(x, altExp_name), "logcounts") <- calc_clr

  return(x)
}
