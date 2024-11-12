#' check_fastq_type
#'
#' @description Determine scRNA-seq fastq type. Fastq can be one of: UMI+Barcode, Sample Index, or Reads. Only works with short-read sequencing.
#'
#'
#' @param x full path to fastq file
#' @param chemistry in development. specify cellRanger chemistry type, v1 - v5, if known for more precise detection.
#' @param n number of rows of fastq file to use
#' @return character string. One of barcode, sample, or read.
#' @import dplyr
#' @import purrr
#' @export

check_fastq_type <- function(x,
                             chemistry = NULL,
                             n = 50){
  if(n < 10){
    warning("n less than 6 is not recommended as there may be limited information to determine single-cell sequencing fastq type")
  }

  # import fastq
  fastq <- readLines(x, n = n)

  # sequence analysis
    # find lines that contain only alpha characters
    # we are ok with quality fields too because length is also representative of sequence length
    # then calculate average length
  seq <- fastq %>% keep(~str_detect(., "^[:alpha:]*$"))
  aveSeqLen <- seq %>% lapply(nchar) %>% unlist %>% mean

  # assign type
  type <- case_when(8 <=aveSeqLen & aveSeqLen <= 10  ~ "sample",
                    26 <=aveSeqLen & aveSeqLen<= 30 ~ "barcode",
                    50 <= aveSeqLen       ~ "read")

  if(is.na(type)) stop("fastq type could not be determined. Sequence length does not match expected sizes for barcode + UMI (26 - 30bp),  sample (8 - 10bp), or read (>=50")

  return(type)




}


