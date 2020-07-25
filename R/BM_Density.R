#' Beta Mixture pdf
#'
#' returns the density of a given point based on the given parameters of a beta mixture distribution
#' @param website the url of a player's page on eliteprospects.com.
#' @return data frame (or list of data frames) with data from all the players drafted in the given year.
#' @export

BM_Density <- function(x, Mixture_Parameters, Alphas, Betas, Log = F) {
  # TODO: add pre-condition checks for Mixture_Parameters, Alphas, and Betas
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
