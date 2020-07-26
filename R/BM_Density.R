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
  #checking x ----
  if (any(typeof(x) != "double")) {stop("x must contain only numeric elements")}
  if (any(x < 0 | x > 1)) {stop("x not between 0 and 1")}
  #checking Alphas and Betas ----
  len = length(Alphas)
  if (len != length(Betas)) {stop("Alphas and Betas must be the same length")}
  if (any(typeof(Alphas) != "double")) {stop("Alphas must contain only numeric elements")}
  if (any(typeof(Betas) != "double")) {stop("Betas must contain only numeric elements")}
  #checking Mixture_Parameters ----
  if (any(typeof(Mixture_Parameters) != "double")) {stop("Mixture_Parameter must contain only numeric elements")}
  if (length(Mixture_Parameters) == len) {
    if (abs(sum(Mixture_Parameters) - 1) > 0.001) {stop("Fully defined Mixture_Parameter must sum to 1")}
  }
  else if (length(Mixture_Parameters) == len-1) {
    if (sum(Mixture_Parameters) > 1 | sum(Mixture_Parameters) < 0) {stop("Partially defined Mixture_Parameter must sum to between 0 and 1")}
  }
  else {stop("Mixture_Parameter must be either the same length or one smaller compared to Alphas/Betas")}
  #checking log ----
  if (typeof(Log) != 'logical') {stop("Log must be a boolean")}



  #actual code to compute mixture density for given points ----
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
