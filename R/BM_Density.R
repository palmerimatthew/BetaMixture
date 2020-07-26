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
  if (any(Mixture_Parameters < 0 | Mixture_Parameters > 1)) {stop("Mixture_Parameters must be between 0 and 1")}
  if (length(Mixture_Parameters) == len) {
    if (abs(sum(Mixture_Parameters) - 1) > 0.001) {stop("Fully defined Mixture_Parameter must sum to 1")}
  }
  else if (length(Mixture_Parameters) == len-1) {
    if (sum(Mixture_Parameters) > 1 | sum(Mixture_Parameters) < 0) {stop("Partially defined Mixture_Parameter must sum to between 0 and 1")}
    temp = sum(Mixture_Parameters)
    Mixture_Parameters = c(Mixture_Parameters, 1-temp)
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



test_error <- function(evaluate_statement, expected_result) {
  temp = evaluate(evaluate_statement)[[2]]
  if (length(attr(temp, "class")) == 0) {
    print(paste("Unexpected Result in", evaluate_statement))
    return(FALSE)
  }
  if (attr(temp, "class")[1] == 'simpleError') {
    if(temp$message != expected_result) {
      print(paste("Expected:", expected_result))
      print(paste("Actual:", temp$message))
      return(FALSE)
    } else {return(TRUE)}
  }
}

Within_error_bounds <- function(BM_dens_value, expression_evaluation, threshold) {
  return(all(abs(BM_dens_value - expression_evaluation) < threshold))
}

#requires package "evaluate"
Test_BM_Density <- function() {
  #testing x checks ----
  x_check1 = test_error("BM_Density(c('a', 1, 2, 3), 1, 5, 6)",
                        "x must contain only numeric elements") #numeric elements
  x_check2 = test_error("BM_Density(c(rep(1, 10000), 'a'), 1, 5, 6)",
                        "x must contain only numeric elements") #numeric elements
  x_check3 = test_error("BM_Density(c(3), 1, 5, 6)",
                        "x not between 0 and 1") #between 0 and 1
  x_check4 = test_error("BM_Density(c(-1, 2), 1, 5, 6)",
                        "x not between 0 and 1") #between 0 and 1
  if(all(x_check1, x_check2, x_check3, x_check4)) {print("Passed all 'x' Tests")}

  #testing Alphas and Betas ----
  AB_test1 = test_error("BM_Density(0.5, 1, c(1,2), 6)",
                        "Alphas and Betas must be the same length") #Alphas longer than Betas
  AB_test2 = test_error("BM_Density(0.5, 1, 6, c(1,2))",
                        "Alphas and Betas must be the same length") #Betas longer than Alphas
  AB_test3 = test_error("BM_Density(0.5, 1, c(rep(1, 10000), 'a'), rep(1, 10001))",
                        "Alphas must contain only numeric elements") #numeric elements
  AB_test4 = test_error("BM_Density(0.5, 1, rep(1, 10001), c(rep(1, 10000), 'a'))",
                        "Betas must contain only numeric elements") #numeric elements
  if(all(AB_test1, AB_test2, AB_test3, AB_test4)) {print("Passed all 'AB' Tests")}

  #testing Mixture_Parameters ----
  MP_test1 = test_error("BM_Density(0.5, c(0.5, 'a'), 1, 2)",
                        "Mixture_Parameter must contain only numeric elements") #numeric elements
  MP_test2 = test_error("BM_Density(0.5, c(-1, 0.5), 1, 2)",
                        "Mixture_Parameters must be between 0 and 1") # < 0
  MP_test3 = test_error("BM_Density(0.5, c(0.5, 2), 1, 2)",
                        "Mixture_Parameters must be between 0 and 1") # > 1
  MP_test4 = test_error("BM_Density(0.5, c(0.5, 0.25), c(1, 2), c(2, 1))",
                        "Fully defined Mixture_Parameter must sum to 1") #fully defined < 1
  MP_test5 = test_error("BM_Density(0.5, c(0.5, 0.5, 0.5), c(1,1,1), c(2,2,2))",
                        "Fully defined Mixture_Parameter must sum to 1") #full defined > 1
  MP_test6 = test_error("BM_Density(0.5, c(0.5, 0.75), c(1,1,1), c(1,1,1))",
                        "Partially defined Mixture_Parameter must sum to between 0 and 1") #partially defined between 0 and 1
  if(all(MP_test1, MP_test2, MP_test3, MP_test4, MP_test5, MP_test6)) {print("Passed all 'MP' Tests")}

  #testing Log ----
  Log_test1 = test_error("BM_Density(0.5, c(0.5, 0.49999), c(1,1), c(1,1), 'a')",
                         "Log must be a boolean") #fully defined Mixture_Parameter sum close to 1 (want to make sure it gets past MP tests, and fails Log test)
  Log_test2 = test_error("BM_Density(0.5, c(0.5, 0.50001), c(1,1), c(1,1), 'b')",
                         "Log must be a boolean") #fully defined Mixture_Parameter sum close to 1 (want to make sure it gets past MP tests, and fails Log test)
  Log_test3 = test_error("BM_Density(0.5, c(0.5, 0.25), c(1,1,1), c(2,2,2), 'a')",
                         "Log must be a boolean") #partially defined Mixture_Parameter between 0 and 1 (want to make sure it gets past MP tests, and fails Log test)
  if(all(Log_test1, Log_test2, Log_test3)) {print("Passed all 'Log' Tests")}

  #testing valid parameters ----
  #fully defined Mixture_Parameters
  Valid_test1 = Within_error_bounds(BM_Density(0.5, c(0.4, 0.6), c(2, 10), c(5, 1), Log = F),
                                    0.4*dbeta(0.5, 2, 5) + 0.6*dbeta(0.5, 10, 1), 0.001)
  Valid_test2 = Within_error_bounds(BM_Density(0.5, c(0.4, 0.6), c(2, 10), c(5, 1), Log = T),
                                    log(0.4*dbeta(0.5, 2, 5) + 0.6*dbeta(0.5, 10, 1)), 0.001)
  #Partially defined Mixture_Parameters
  Valid_test3 = Within_error_bounds(BM_Density(c(0.25, 0.75, 0.5), c(0.3, 0.6), c(1, 2, 3), c(10, 5, 15), Log = F),
                                    0.3*dbeta(c(0.25, 0.75, 0.5), 1, 10) + 0.6*dbeta(c(0.25, 0.75, 0.5), 2, 5) + 0.1*dbeta(c(0.25, 0.75, 0.5), 3, 15), 0.001)

}
