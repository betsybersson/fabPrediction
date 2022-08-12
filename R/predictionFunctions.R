#' Obtain a FAB conformal prediction interval
#'
#' This function computes a FAB conformal prediction region as described in 
#' Bersson and Hoff 2022.
#'
#' @param Y Observed data vector
#' @param mu Prior expected mean of the population mean 
#' @param tau2 Prior expected variance of the population mean
#' @param alpha Prediction error rate
#' @return alpha-level prediction interval
#' @export
fab_conf_pred = function(Y,mu = 0,tau2 = 1,alpha = .05){
  
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
  
  # if ( k == 0 ) {
  #   stop("ALPHA is set too low.")
  # }
  
  # return interval
  int # return(list("bounds" = int, "prob" = (1-k/(N+1))*100))
  
}

#' Obtain a distance-from-average conformal prediction interval
#'
#' This function computes a conformal prediction region under the distance-from-average
#' non-conformity measure. That is, |a + bz*| <= |ci + di z^*| where i indexes training data.
#'
#' @param Y Observed data vector
#' @param alpha Prediction error rate
#' @return alpha-level prediction interval
#' @export
avg_conf_pred = function(Y,alpha){

  ## get helpers
  N = length(Y)
  sumY = sum(Y)
  avgY = (N+1) * Y
  
  a = sumY
  b = -N
  ci = sumY-avgY # constant
  di = rep(1,N) # multiplied by zstar
  
  b1 = (a-ci)/(di-b)
  b2 = (-a-ci)/(di+b)
  
  # obtain final solution
  S = sort(c(b1,b2))
  k = floor(alpha*(N+1))
  int = c(S[k],S[2*N-k+1])
  
  # if ( k == 0 ) {
  #   stop("ALPHA is set too low.")
  # }
  
  # return interval
  int #return(list("bounds" = int, "prob" = (1-k/(N+1))*100))
  
}
