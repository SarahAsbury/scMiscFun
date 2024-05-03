#' split_sce
#'
#' split sce object into multiple objects by variable
#'
#' @param sce sce object
#' @param by character string; variable to split by
#' @return list of sce objects
#' @examples
#' split_sce(sce = sce, by = "sample_id")
#' @import dplyr
#' @import purrr
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @export

split_sce <- function(sce,
                      by
)
{

  group <- sce %>% pull(by) %>% unique

  map(group,

      ~sce %>% dplyr::filter(!!sym(by) == .x)
  ) %>%
    setNames(group)

}
