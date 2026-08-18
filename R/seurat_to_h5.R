#' seurat_to_h5
#'
#' Converts RNA and PROTEIN experiments from a seurat object into an h5ad file for import to Python.
#'
#' @param x seurt object
#' @param output_dir_prefix # full path to output directory + any prefix. Will append assays and h5.Seurat to the file
#' @param assays sce object assay name to  to export. one file will be exported per each experiment.
#' @param layer optional layer name to export to X in H5AD. If NULL (default), exports counts.
#' @return exports h5seurat and h5ad files from a seurat object
#' @examples seurat_to_h5(seurat, "path/to/output/directory/output_prefix")
#' @examples seurat_to_h5(seurat, "path/to/output/directory/output_prefix", layer = "tcell_denoised_protein")
#' @import dplyr
#' @import purrr
#' @import Seurat
#' @import SeuratDisk
#' @export


seurat_to_h5 <- function(x, #seurat
                         output_dir_prefix,
                         assays = c("RNA", "PROTEIN"),
                         layer = NULL  # if NULL, exports counts
){

  map(assays,

      \(assay){

        # If a specific layer is requested, check if it exists before proceeding
        if (!is.null(layer) && layer != "counts") {
          if (!layer %in% Assays(x)) {
            cat(sprintf("  Layer '%s' not found in Seurat object, skipping %s export\n", layer, assay))
            return(NULL)
          }
        }

        # params - add layer to filename if specified
        layer_suffix <- if (!is.null(layer) && layer != "counts") paste0("_", layer) else ""
        output_fn <- sprintf("%s_%s%s.h5Seurat", output_dir_prefix, assay, layer_suffix)
        cat("\n export h5ad file:", output_fn, "\n\n")

        # dm
        export_seurat <- x
        DefaultAssay(export_seurat) <- assay

        # If a different layer is specified, swap it into counts slot
        if (!is.null(layer) && layer != "counts") {
          cat(sprintf("  Exporting layer '%s' to X in %s H5AD\n", layer, assay))
          export_seurat[[assay]]@counts <- export_seurat[[layer]]@counts
        }

        export_seurat[[assay]]$data <- NULL  # remove log counts to export counts

        # export h5Seurat
        SaveH5Seurat(export_seurat,
                     filename = output_fn,
                     overwrite = T)

        # convert to h5ad
        Convert(output_fn,
                dest = "h5ad",
                overwrite = T)
      }
  )

}


