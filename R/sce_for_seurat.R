#' sce_for_seurat
#'
#' converts RNA and PROTEIN assays from single cell experiment to seurat format. reduced dimensions are removed by default unless specified to be kept.
#'
#' @param x sce object
#' @param keep_dim characteor vector of dimensions to keep
#' @return seurat object
#' @examples sce_for_seurat(keep_dim = "tcell_totalVI")
#' @import dplyr
#' @import purrr
#' @import Seurat
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @export



sce_for_seurat <- function(x,
                           keep_dim = NULL)
{

  # clean
  # remove extra assays
  assay(x, "rel_count") <- NULL

  # remove extra altExp
  remove_exp <- altExpNames(x) %>% str_subset("^RNA$|^PROTEIN$", negate = T)
  for(i in remove_exp){
    altExp(x, i) <- NULL
  }

  # remove reduced dim
  for(i in reducedDimNames(x)[!(reducedDimNames(x) %in% keep_dim)]){
    reducedDim(x, i) <- NULL
  }

  # format raw counts for seurat
  assay(x, "counts") <-  assay(x, "rawcounts") %>% as("sparseMatrix")
  assay(x, "logcounts") %<>% as("sparseMatrix")

  counts(altExp(x, "PROTEIN")) <- altExp(x, "PROTEIN") %>% assay("rawcounts") %>% as("sparseMatrix")
  logcounts(altExp(x, "PROTEIN")) %<>% as("sparseMatrix")

  return(x)
}
