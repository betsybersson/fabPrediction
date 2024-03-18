#' Obtain a distance-to-average conformal prediction interval
#'
#' This function computes a conformal prediction region under the distance-from-average
#' non-conformity measure. That is, |a + bz*| <= |ci + di z^*| where i indexes training data.
#'
#' @param Y Observed data vector
#' @param alpha Prediction error rate
#' @return  pred object
#' @imports package
#' @export
dtaPrediction = function(Y,alpha = .15){
  
  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    print("Y converted to vector!")
  }
  
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
  
  out = list("bounds" = int, "coverage" = (1-k/(N+1))*100,
             "data" = Y, "class" = "continuous")
  
  class(out) = 'pred'
  
  # return pred object
  return(out)
}