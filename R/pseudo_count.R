#' Add pseudocounts to numeric data
#'
#' Adds half of the smallest non-zero value in the dataset to each value. This is
#' useful for handling zero values before log-ratio transformations.
#'
#' @param x Numeric vector containing the data
#'
#' @return Numeric vector with pseudocounts added
#'
#' @examples
#' data <- c(0, 1, 2, 3, 0, 4)
#' pseudo_count(data)
#'
#' @export
pseudo_count <- function(x) {
  x + min(x[x != 0])/2
}
