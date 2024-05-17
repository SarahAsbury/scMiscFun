#' plot_sample_qc
#'
#' plot standard quality control metrics of total UMI and total detected genes per sample
#'
#' @param sce sce object
#' @param plot_by character vector; variables to plot by
#' @return list of plots
#' @examples
#' plot_sample_qc(sce, plot_by = c("sample", "annot"))
#' @import dplyr
#' @import purrr
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @import scater
#' @import cowplot
#' @export

plot_sample_qc <- function(sce,
                           plot_by){
  map(c(plot_by),

      ~plot_grid(

        # sequencing depth + detected genes
        plotColData(sce, x = "sum", y="detected", colour_by=.x) +
          theme(legend.position = "none") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Total UMI",
               y = "Total genes"),

        plot_grid(plotColData(sce, x = .x, y="sum", colour_by=.x) +
                    theme(legend.position = "none") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs(y = "Total UMI"),

                  plotColData(sce, x = .x, y="detected", colour_by=.x) +
                    theme(legend.position = "none") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs(x = "Total genes"),

                  ncol = 2
        ),


        ncol = 1
      )
  ) %>% setNames(plot_by)

}
