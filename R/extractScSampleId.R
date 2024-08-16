#' extractScSampleId
#'
#' extract sample Id from common single cell formats
#'
#' @param x filename to extract sample Id from
#' @param auto binary. if true, will attempt to extract common scrnaseq sample Id formats from filename
#' @param type character vector. specify which common sample IDs from common single cell formats to extract. Options include: sraSample, sraRun, sraExperiment, sraProject, sraAccession
#' @param regex user specified regex to extract sampleId. If provided alongside type, any match to common regex type(s) or user-input regex will be used
#' @return sample ID
#' @import dplyr
#' @import purrr
#' @export

extractScSampleId <- function(x,
                              auto = F,
                              type = NULL,
                              regex = NULL){

  # --- regex repo ---
  common_regex_patterns <- list(
    sraSample = "SAMN[:digit:]{8,9}|SRS[:digit:]{6,7}",
    sraRun = "SRR[:digit:]{6,9}",
    sraExperiment = "SRX[:digit:]{6,9}",
    sraProject = "SRP[:digit:]{6,8}",
    sraAccession = "SRA[:digit:]{6,8}"
  )

  # --- generate working regex ---
  working_regex_pre <- list()
  if(auto){
    # autodetect common regex type
    auto_regex = map(common_regex_patterns, ~str_detect(x, .x)) %>% .[. == T] %>% names
    if(length(auto_regex) > 1){
      stop("Auto detection failed because more than 1 common regex pattern was detected in the filename. Please manually specify the sample name or preprocess the filename so only one type can be included at a time.")
    }

    working_regex_pre <- c(working_regex_pre, common_regex_patterns[names(common_regex_patterns) %in% auto_regex])
  }

  if(!is.null(type)){
    # user-specified common regex type
    working_regex_pre <- c(working_regex_pre, common_regex_patterns[names(common_regex_patterns) %in% type])
  }

  if(!is.null(regex)){
    # user-specified regex
    working_regex_pre <- c(working_regex_pre, list(user_regex = regex))
  }

  working_regex <- working_regex_pre %>% paste(collapse = "|")


  # --- extract ---
  out <- str_extract(x, working_regex)



  return(out)

}

# tests
# filename <- c("/ddn_exa/campbell/sasbury/sasbury/fastq-dump/data/lung/icb/cellranger/count/SAMN35370102/emptydrops_filtered")
# extractScSampleId(filename, type = "sraSample")
# extractScSampleId(filename, auto = T)
# extractScSampleId(filename, regex = "SAMN[:digit:]{8,9}")



