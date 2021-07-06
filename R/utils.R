#' Compute z-score
#'
#' @param x A numeric vector
#' @param na.rm Should NA values be excluded [default: TRUE]
#' @return A vector as long as \code{x} with the mean and standard deviation
#'   normalized.
Z <- function(x, na.rm = TRUE) {
    return((x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm))
}

#' Compute standard error of mean
#'
#' @param x A numeric vector
#' @param na.rm Should NA values be excluded [default: FALSE]
#' @return A scalar value corresponding to the standard error of the mean of x.
se <- function(x, na.rm = FALSE) {
    if (na.rm == TRUE) {
        z <- !is.na(x)
        x <- x[z]
        n <- sum(z)
    } else {
        n <- length(x)
    }
    return(sd(x) / sqrt(n))
}

#' Convert 1D indexes into 2D subscript (i.e., coordinate)
#' 
#' @param ind A vector of indexes to convert
#' @param m The number of rows in the matrix the subscript should respect.
#' @param for_igraph A logical. If true, a vector is returned, such that each
#'   consecutive pair of elements represent a 2D coordinate. [default: FALSE]
#' @return Coordinate pairs either as a matrix with two columns or as a vector
#'   (depending on \code{for_igraph}).
ind2sub <- function(ind, m, for_igraph = FALSE) {
  r = ((ind - 1) %% m) + 1
  c = floor((ind-1) / m) + 1
  x <- if (for_igraph) c(rbind(r, c)) else cbind(r, c)
  return(x)
}
