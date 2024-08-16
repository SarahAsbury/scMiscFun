#' get_symbol
#'
#' @description query a vector of gene ids using biomaRt ensembl database and return symbol. Optionally return additional attributes.
#'
#' @param gene_ids character vector of gene IDs (required)
#' @param from_type Optional. specify format of gene IDs. One of:  If not provided, will automatically try and detect format.
#' @param add_attribute Optional. additional gene attributes to add. use useMart("ensembl", dataset = "hsapiens_gene_ensembl") %>% listAttributes to see what attributes are available.
#' @return dataframe of gene symbol and any additional attributes
#' @import dplyr
#' @import purrr
#' @import biomaRt
#' @export



get_symbol <- function(gene_ids, from_type = "auto", add_attribute = NULL) {
  # use ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  # auto detect
  if(from_type == "auto"){

    from_type <- detect_gene_format(gene_ids)
  }

  # Map the input type to the correct BioMart attribute
  input_attr <- switch(from_type,
                       "ensembl" = "ensembl_gene_id",
                       "entrez" = "entrezgene_id",
                       "refseq" = "refseq_mrna",
                       "hgnc" = "hgnc_id",
                       "ucsc" = "ucsc",
                       stop("Unsupported input type"))

  attributes <- c(input_attr, "external_gene_name", add_attribute)

  # Perform the query
  result <- getBM(attributes = attributes,
                  filters = input_attr,
                  values = gene_ids,
                  mart = ensembl
  ) %>%
    mutate(
      Symbol = ifelse(external_gene_name != "",
                      external_gene_name,
                      !!sym(input_attr)
      )
    ) %>%
    column_to_rownames(var = input_attr) %>%
    dplyr::select(Symbol, any_of(add_attribute))

  return(result)
}
