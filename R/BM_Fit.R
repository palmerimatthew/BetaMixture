#' Beta Mixture Density Fitting Procedure
#'
#' Fits a Beta Mixture density to a given set of data
#' @param x vector of data points. Will fit beta mixture density to this.
#' @param K Number of distributions composing mixture
#' @param threshold Determines "accuracy" of final returned parameters
#' @param seed Sets seed for random processes. Should only be used if replicable results are needed.
#' @return TODO
#' @export
BM_Fit <- function(x, K, threshold, seed = NA) {
  #setting seed if seed is not NA
  if (!is.na(seed)){set.seed(seed)}


}
