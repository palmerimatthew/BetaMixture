#' Beta Mixture pdf
#'
#' returns the density of a given point based on the given parameters of a beta mixture distribution
#' @param x vector of data points. Will compute the density of them.
#' @param Mixture_Parameters Vector of mixture parameters. Must add up to 1 if of length N, or must be within (0, 1) if of length N-1 (reference Alphas and Betas for help with length).
#' @param Alphas Vector of Alpha shape parameters. Must be of length N (reference Alphas and Betas for help with length).
#' @param Betas Vector of Beta shape parameters. Must be of length N (reference Alphas and Betas for help with length).
#' @param Log Boolean to deterime if the log of the density is to be returned.
#' @return Vector of density measures. Same length as input vector x, with cooresponding indices.
#' @export

BM_Density <- function(x, Mixture_Parameters, Alphas, Betas, Log = F) {
  # TODO: add pre-condition checks for Mixture_Parameters, Alphas, and Betas

  #actual code to compute mixture density for given points
  len = length(Alphas)
  dens = numeric(length(x))
  for (i in 1:len) {
    dens = dens + Mixture_Parameters[i]*dbeta(x, Alphas[i], Betas[i], log = F)
  }
  if (Log) {
    return (log(dens))
  } else {
    return(dens)
  }
}
