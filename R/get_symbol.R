#' get_symbol
#'
#' @description query a vector of gene ids using biomaRt ensembl database and return symbol. Optionally return additional attributes.
#'
#' @param gene_ids character vector of gene IDs (required)
#' @param from_type Optional. specify format of gene IDs. One of:  If not provided, will automatically try and detect format.
#' @param add_attribute Optional. additional gene attributes to add. use listAttributes for ensembl & hsapiens_gene_ensembl to see what attributes are available.
#' @param missing Default = TRUE. Binary. If TRUE, any gene IDs not found in the ensembl database will still be included in the results. Symbol will be the original gene ID instead.
#'
#' @return dataframe of gene symbol and any additional attributes
#'
#' @import dplyr
#' @import purrr
#' @import biomaRt
#'
#' @export



get_symbol <- function(gene_ids, from_type = "auto", add_attribute = NULL, missing = T) {
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

  # add any missing query
  if(missing == T){
    missing_query <- anti_join(
      data.frame(query = gene_ids),
      result %>% rownames_to_column(var = "query"),
      by = "query"
    ) %>%
      mutate(
        Symbol = query
      ) %>%
      column_to_rownames(var = "query")

    result <- dplyr::bind_rows(result, missing_query)
  }


  return(result)
}
