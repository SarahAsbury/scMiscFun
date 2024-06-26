#' sce_for_seurat2
#'
#' converts assays from single cell experiment to seurat format. Designed for h5 export for scvi. Reduced dimensions are removed by default unless specified to be kept.
#'
#' @param x sce object
#' @param rawcounts name of assay containing rawcounts
#' @param keep_additional_assays character vector of any additional assays - other than rawcounts - to include in seurat object
#' @param altExp which alt experiments to include in seurat object
#' @param keep_dim character vector of dimensions to keep
#' @return seurat object. count layer is the rawcounts layer specified.
#' @examples sce_for_seurat(keep_dim = "tcell_totalVI")
#' @import dplyr
#' @import purrr
#' @import Seurat
#' @import stringr
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @export



sce_for_seurat2 <- function(x,
                            rawcounts = "rawcounts",
                            keep_additional_assays = NULL,
                            keep_dim = NULL,
                            altExp = NULL
                            )
{
  # setup
  regex <- list()

  # clean
  # remove extra assays
  regex$assays <- paste0("^", c(rawcounts, keep_additional_assays), "$", collapse = "|")
  remove_assays <- assayNames(x) %>% str_subset(regex$assays, negate = T)
  for(i in remove_assays){
    assay(x, i) <- NULL
  }

  # remove extra altExp
  regex$altExp <- paste0("^", altExp, "$", collapse = "|")
  remove_exp <- altExpNames(x) %>% str_subset(regex$altExp, negate = T)
  for(i in remove_exp){
    altExp(x, i) <- NULL
  }

  # remove reduced dim
  for(i in reducedDimNames(x)[!(reducedDimNames(x) %in% keep_dim)]){
    reducedDim(x, i) <- NULL
  }

  # format
  # format default experiment assays for seurat
  assay(x, "counts") <-  assay(x, rawcounts) %>% as("sparseMatrix")
  for(i in keep_additional_assays){
    assay(x, i) %<>% as("sparseMatrix")
  }

  # format alt experiments
  if(!is.null(altExp)){
    warning("sce_for_seurat2 has not been tested with altExp. Remove this message once testing occurs.")
    for(i in altExp){
      counts(altExp(x, i)) <- altExp(x, i) %>% assay(rawcounts) %>% as("sparseMatrix")
    }
    for(j in keep_additional_assays){
      assay(altExp(x, "PROTEIN"), i) %<>% as("sparseMatrix")
    }


  }

  return(x)
}
