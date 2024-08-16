#' detect_gene_format
#'
#' @description given a vector of gene ids, automatically detect the gene fomrat from standard formats (ensembl, entrez, refseq, hgnc, uscs). Only provides with format declaration if all genes IDs match that format.
#'
#' @param x character vector of gene IDs
#' @return format of gene IDs
#' @import dplyr
#' @import purrr
#' @import biomaRt
#' @export


detect_gene_format <- function(gene_ids){

  huGeneFormatRegex = list(
    ensembl = "^ENSG[:digit:]+$",
    entrez = "^[:digit:]+$",
    refseq = "^NG_[:digit:]+$",
    hgnc = "^HGNC:[:digit:]+$",
    ucsc = "^hg[:digit:]{1,2}\\..*$"
  )

  checkType <- map(huGeneFormatRegex,

                   function(regex){
                     detectRegex <- gene_ids %>% map(~str_detect(.x, regex)) %>% unlist
                     detectRegexAll <- all(detectRegex)
                     detectRegexAny <- any(detectRegex)

                     return(c(all = detectRegexAll, any = detectRegexAny))
                   }
  )

  if(length(checkType %>% keep(~.x[["all"]])) == 1){
    # if one complete match, use it
    geneFormat <- checkType %>% keep(~.x[["all"]]) %>% names

  } else if (length(checkType %>% keep(~.x[["any"]])) > 0){
    # checks for cases of partial match and inform user
    print(checkType)
    stop("Either multiple matches were identified or some gene IDs did not match the auto-detected format")
  } else{
    # fail
    print(checkType)
    stop("Gene IDs did not match standard formats.")
  }

  return(geneFormat)

}



