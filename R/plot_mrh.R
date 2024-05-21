#' plot_sample_qc
#'
#' plot standard quality control metrics of mitochondrial, ribosomal, and hemoglobin gene counts. designed to work with run_sc_qc pipeline only.
#'
#' @param sce sce object
#' @param plot_by character vector; variables to plot by
#' @param mit_cutoff numeric value; a horizontal line will be added at this threshold.
#' @param include character vector containing m, r, and/or h. If m is in vector, mitochondrial genes will be plotted. If r, ribosomal genes. If h, hemoglobin genes.
#' @return list of plots
#' @examples
#' plot_mrh(sce, plot_by = c("sample", "annot"), mit_cutoff = 5, include = "mrh")
#' @import dplyr
#' @import purrr
#' @import SingleCellExperiment
#' @import tidySingleCellExperiment
#' @import scater
#' @import cowplot
#' @export

plot_mrh <- function(sce,
                     plot_by,
                     mit_cutoff,
                     include = "mrh"
                     )
  {
  map(c(plot_by),

      ~{
        plots <- list()
        if(str_detect(include, regex("m", ignore_case = T))){
          plots <- c(plots,
                     list(
                       plotColData(sce, x = "sum", y="subsets_Mito_percent", colour_by=.x) +
                         theme(legend.position = "none") +
                         geom_hline(yintercept = mit_cutoff, linetype="dashed", color = "black", size = 0.4),

                       plotColData(sce, x = .x, y="subsets_Mito_percent", colour_by=.x) +
                         theme(legend.position = "none") +
                         geom_hline(yintercept = mit_cutoff, linetype="dashed", color = "black", size = 0.4) +
                         theme(axis.text.x = element_blank())
                     )
          )
        }

        if(str_detect(include, regex("r", ignore_case = T))){
          plots <- c(plots,
                     list(
                       plotColData(sce, x = "sum", y="subsets_Ribo_percent", colour_by=.x) +
                         theme(legend.position = "none"),

                       plotColData(sce, x = .x, y="subsets_Ribo_percent", colour_by=.x) +
                         theme(legend.position = "none") +
                         theme(axis.text.x = element_blank())
                     )
          )

        }

        if(str_detect(include, regex("h", ignore_case = T))){
          plots <- c(plot,
                     list(plotColData(sce, x = "sum", y="subsets_Hemo_percent", colour_by=.x) +
                            theme(legend.position = "none"),

                          plotColData(sce, x = .x, y="subsets_Hemo_percent", colour_by=.x) +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_blank())
                     )
          )
        }




        plot_grid(plotlist = plots,
                  ncol = 2
        ) %>%
          # add legend to plot
        plot_grid(.,
                  plotColData(sce, x = .x, y="sum", colour_by=.x) %>% cowplot::get_legend(),
                  rel_widths = c(5,1),
                  ncol = 2)
      }) %>%
            setNames(plot_by)
}

