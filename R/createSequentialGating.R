#' Create an Interactive Sequential Gating Application
#'
#' @description
#' Creates a Shiny application for sequential cell gating analysis of single-cell data.
#' The application allows users to iteratively apply gates to their data, visualize results,
#' and export the final gated population and gate coordinates.
#'
#' @param sce_initial A SingleCellExperiment object containing the initial dataset to be gated
#' @param marker_pairs A list of vectors, where each vector contains four elements:
#'        \itemize{
#'          \item First marker name (x-axis)
#'          \item Second marker name (y-axis)
#'          \item Gate signs (determining positive quadrant)
#'          \item Population name
#'        }
#' @param assay_name Character string specifying which assay in the SingleCellExperiment to use
#'
#' @return A Shiny application object that when run provides:
#'   \itemize{
#'     \item Interactive scatter plots for each gating step
#'     \item Adjustable gate positions using sliders
#'     \item Statistics on cells passing each gate
#'     \item Options to download:
#'       \itemize{
#'         \item Gates as CSV file ("gates.csv")
#'         \item Gated SingleCellExperiment object ("sce_gated.rds")
#'         \item Gate visualizations as images
#'       }
#'   }
#'
#' @details
#' The function creates a sequential gating interface where each gate becomes available
#' only after the previous gate has been applied. For each gate, users can:
#' \itemize{
#'   \item Adjust vertical and horizontal gate positions using sliders
#'   \item View the scatter plot of cells with current gate positions
#'   \item See statistics on the percentage of cells passing the gate
#'   \item Apply the gate and move to the next population
#' }
#'
#' The final gated population can be downloaded along with gate coordinates and visualizations.
#'
#' @examples
#' \dontrun{
#' # Define marker pairs for gating
#' marker_pairs <- list(
#'   c("CD3", "CD19", "++", "T cells"),
#'   c("CD4", "CD8", "+-", "CD4 T cells")
#' )
#'
#' # Create and run the gating application
#' createSequentialGating(
#'   sce_initial = your_sce_object,
#'   marker_pairs = marker_pairs,
#'   assay_name = "protein"
#' )
#' }
#'
#' @import shiny
#' @import ggplot2
#' @import SingleCellExperiment
#' @import magrittr
#' @import dplyr
#' @import purrr
#' @import tidySingleCellExperiment
#'
#' @export


createSequentialGating <- function(sce_initial, marker_pairs, assay_name) {
  ui <- fluidPage(
    imap(
      marker_pairs,
      function(markers, i) {
        panel_name <- paste0("Gate ", i)

        conditionalPanel(
          condition = if(i == 1) "true" else paste0("input.apply_gate_", i-1, " > 0"),
          h3(panel_name),
          fluidRow(
            column(6,
                   sliderInput(paste0("vline_", i), "Vertical Gate Position:",
                               min = min(assay(sce_initial, assay_name)[markers[1], ]),
                               max = max(assay(sce_initial, assay_name)[markers[1], ]),
                               value = 0, step = 0.1)
            ),
            column(6,
                   sliderInput(paste0("hline_", i), "Horizontal Gate Position:",
                               min = min(assay(sce_initial, assay_name)[markers[2], ]),
                               max = max(assay(sce_initial, assay_name)[markers[2], ]),
                               value = 0, step = 0.1)
            )
          ),
          plotOutput(paste0("plot_", i), height = "600px"),
          verbatimTextOutput(paste0("stats_", i)),
          # Add Done button for all except last gate
          if(i < length(marker_pairs)) {
            actionButton(paste0("apply_gate_", i), "Apply Gate and Continue")
          } else {
            # Add download gates when done
            tagList(
              actionButton("download_gates", "Download Gates"),
              actionButton("exit_app", "Exit")
            )
          }
        )
      })
  )

  server <- function(input, output, session) {
    # initialize reactive variables
    ## filtered_sces
    filtered_sces <- reactiveValues()
    filtered_sces$data <- list(sce_initial)
    ## boundaries
    boundaries <- reactiveValues()
    boundaries$data <- list()

    # Create plots and stats
    iwalk(
      marker_pairs,
      function(markers, i) {

        # iterative plotting
        output[[paste0("plot_", i)]] <- renderPlot({

          # get current/filtered sce
          current_sce <- filtered_sces$data[[i]]

          # plot
          create_gate_plot(
            current_sce,
            marker1 = markers[1],
            marker2 = markers[2],
            x_boundary = input[[paste0("vline_", i)]],
            y_boundary = input[[paste0("vline_", i)]],
            population_name =  markers[4],
            assay_name = assay_name)
        })


        #   plot_data <- data.frame(
        #     x = assay(current_sce, assay_name)[markers[1], ],
        #     y = assay(current_sce, assay_name)[markers[2], ]
        #   )
        #
        #   # plot
        #   ggplot(plot_data, aes(x = x, y = y)) +
        #     geom_point(alpha = 0.2) +
        #     geom_vline(xintercept = input[[paste0("vline_", i)]],
        #                color = "red", linetype = "dashed") +
        #     geom_hline(yintercept = input[[paste0("hline_", i)]],
        #                color = "red", linetype = "dashed") +
        #     theme_minimal() +
        #     labs(title = markers[4],
        #          x = markers[1],
        #          y = markers[2]
        #     )
        # }

        # add statistics
        output[[paste0("stats_", i)]] <- renderText({
          current_sce <- filtered_sces$data[[i]]

          # total number of cells that are in the positive gate (top right quadrant)
          total <- ncol(current_sce)
          passing <- apply_gate_signs(
            sce = current_sce,
            assay_name = assay_name,
            marker1 = markers[1],
            marker2 = markers[2],
            gate_signs = markers[3],
            boundaries = list(
              x = input[[paste0("vline_", i)]],
              y = input[[paste0("hline_", i)]]
            )
          ) %>% sum

          paste0("Positive cells: ", passing, " (", round(100*passing/total, 1), "%)")
        })
      })

    # Handle gate applications
    iwalk(
      marker_pairs,
      function(markers, i) {
        # Safely handle next markers, using NULL for the last gate
        next_markers <- if(i < length(marker_pairs)) marker_pairs[[i + 1]] else NULL

        observeEvent(input[[paste0("apply_gate_", i)]], {
          # get sce
          current_sce <- filtered_sces$data[[i]]

          # get gate boundaries
          boundaries$data[[i]] <- list(
            population = markers[4],
            x = input[[paste0("vline_", i)]],
            y = input[[paste0("hline_", i)]]
          )

          # filter cells
          cells_pass <- apply_gate_signs(
            sce = current_sce,
            assay_name = assay_name,
            marker1 = markers[1],
            marker2 = markers[2],
            gate_signs = markers[3],
            boundaries = boundaries$data[[i]]
          )

          # Only update filtered_sces and sliders if not the last gate
          if(i < length(marker_pairs)) {
            filtered_sces$data[[i + 1]] <- current_sce[, cells_pass]

            # get next sce and set interactive sliders
            next_sce <- filtered_sces$data[[i + 1]]
            updateSliderInput(session, paste0("vline_", i + 1),
                              min = min(assay(next_sce, assay_name)[next_markers[1], ]),
                              max = max(assay(next_sce, assay_name)[next_markers[1], ]))
            updateSliderInput(session, paste0("hline_", i + 1),
                              min = min(assay(next_sce, assay_name)[next_markers[2], ]),
                              max = max(assay(next_sce, assay_name)[next_markers[2], ]))
          }
        })
      })

    # Download results
    observeEvent(input$download_gates, {

      ### save filtered sce
      # filter on last gate
      # get last sce
      last_sce <- filtered_sces$data[[length(filtered_sces$data)]]

      # get last marker pairs
      last_markers <- marker_pairs[[length(marker_pairs)]]

      # get last boundaries
      last_boundaries <- list(
        population = last_markers[4],
        x = input[[paste0("vline_", length(marker_pairs))]],
        y = input[[paste0("hline_", length(marker_pairs))]]
      )

      # filter cells
      cells_pass <- apply_gate_signs(
        sce = last_sce,
        assay_name = assay_name,
        marker1 = last_markers[1],
        marker2 = last_markers[2],
        gate_signs = last_markers[3],
        boundaries = last_boundaries
      )
      output_sce <- last_sce[, cells_pass]

      # download
      if(swap_exp == TRUE){
        output_sce %<>% swapAltExp("RNA")
      }
      saveRDS(
        output_sce,
        file = "sce_gated.rds"
      )

      ### save gates
      gate_df <- map_dfr(
        c(boundaries$data, last_boundaries),

        ~.x %>% data.frame %>% setNames(c("population", "x_gate", "y_gate"))
      )

      write.csv(
        gate_df,
        "gates.csv",
        row.names = FALSE
      )


      ### save report
      screenshot(
        filename = "gates",
        scale = 3,
        download = FALSE,
        server_dir="."
      )

    }
    )

    # exit
    observeEvent(input$exit_app, {
      stopApp()
    }
    )
  }

  shinyApp(ui = ui, server = server)
}
