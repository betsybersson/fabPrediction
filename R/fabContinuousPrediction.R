#' Obtain a FAB conformal prediction interval
#'
#' This function computes a FAB conformal prediction region as described in
#' Bersson and Hoff 2022.
#'
#' @param Y Observed data vector
#' @param mu Prior expected mean of the population mean
#' @param tau2 Prior expected variance of the population mean
#' @param alpha Prediction error rate
#' @return pred object
#' @export
fabContinuousPrediction = function(Y,alpha = .15,mu = 0,tau2 = 1){
  
  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    message("Y converted to vector!")
  }
  
  # parameters
  N = length(Y)
  tau2_theta = 1/(1/tau2 + N + 1)
  
  # critical values
  sol1s = Y
  sol2s = (2*(mu/tau2 + sum(Y))*tau2_theta-Y)/
    (1-2*tau2_theta)
  
  # obtain alpha level prediction interval
  S = sort(c(sol1s,sol2s))
  k = floor(alpha*(N+1))
  int = c(S[k],S[2*N-k+1])
  
  out = list("bounds" = int, "coverage" = (1-k/(N+1))*100,
             "data" = Y, "class" = "continuous")
  
  class(out) = 'pred'
  
  # return pred object
  return(out)
}