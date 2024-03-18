#' Obtain a Bayesian prediction interval
#'
#' This function computes a Bayesian prediction interval based on a normal model.
#'
#' @param Y Observed data vector
#' @param mu Prior expected mean of the population mean
#' @param tau2 Prior expected variance of the population mean
#' @param alpha Prediction error rate
#' @return  pred object
#' @export
bayesNormalPrediction = function(Y,alpha = .15,mu = 0,tau2 = 1){
  
  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    print("Y converted to vector!")
  }
  
  N = length(Y)
  ybar = mean(Y)
  sig2 = var(Y) # use county specific sd
  
  varj = 1 / (1/tau2 + N/sig2)
  thetaj = (mu/tau2 + ybar * N/sig2) * varj
  
  tval = qt(1-alpha/2,df = N-1)
  
  ci_help = tval * sqrt(varj + sig2)
  
  out = list("bounds" = c(thetaj - ci_help,thetaj + ci_help),
             "coverage" = (1-alpha)*100,
             "data" = Y, "class" = "continuous")
  
  class(out) = 'pred'
  
  # return pred object
  return(out)
}