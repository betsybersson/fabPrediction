#' Wrapper to obtain a prediction interval for continuous data
#'
#' This function computes a prediction interval from a number of methods.
#'
#' @param Y Observed data vector
#' @param method Choice of prediction method. Options include FAB, DTA, direct, Bayes.
#' @param epsilon Prediction error rate
#' @param mu Prior expected mean of the population mean
#' @param tau2 Prior expected variance of the population mean
#' @return pred object containing prediction interval bounds and interval coverage
#' @export
predictionInterval = function(Y,method = "FAB",
                          epsilon = .15,
                          mu = 0,tau2 = 1){

  if (method == "FAB"){
    out = fabPrediction(Y,epsilon = epsilon,mu = mu,tau2 = tau2)
  } else if (method == "DTA"){
    out = dtaPrediction(Y,epsilon = epsilon)
  } else if (method == "direct"){
    out = normalPrediction(Y,epsilon = epsilon)
  } else if (method == "Bayes"){
    out = bayesNormalPrediction(Y,epsilon = epsilon,mu = mu,tau2 = tau2)
  } else {
    stop(paste0("Error! Method ",method," is not a valid option!"))
  }

  # return pred object
  return(out)
}
#' Obtain a FAB conformal prediction interval
#'
#' This function computes a FAB conformal prediction region as described in
#' Bersson and Hoff 2022.
#'
#' @param Y Observed data vector
#' @param mu Prior expected mean of the population mean
#' @param tau2 Prior expected variance of the population mean
#' @param epsilon Prediction error rate
#' @return pred object
#' @export
fabPrediction = function(Y,epsilon = .15,mu = 0,tau2 = 1){

  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    print("Y converted to vector!")
  }

  # parameters
  N = length(Y)
  tau2_theta = 1/(1/tau2 + N + 1)

  # critical values
  sol1s = Y
  sol2s = (2*(mu/tau2 + sum(Y))*tau2_theta-Y)/
    (1-2*tau2_theta)

  # obtain epsilon level prediction interval
  S = sort(c(sol1s,sol2s))
  k = floor(epsilon*(N+1))
  int = c(S[k],S[2*N-k+1])

  # if ( k == 0 ) {
  #   stop("ALPHA is set too low.")
  # }

  out = list("bounds" = int, "coverage" = (1-k/(N+1))*100,
             "data" = Y, "class" = "continuous")

  class(out) = 'pred'

  # return pred object
  return(out)
}
#' Obtain a Bayesian prediction interval
#'
#' This function computes a Bayesian prediction interval based on a normal model.
#'
#' @param Y Observed data vector
#' @param mu Prior expected mean of the population mean
#' @param tau2 Prior expected variance of the population mean
#' @param epsilon Prediction error rate
#' @return  pred object
#' @export
bayesNormalPrediction = function(Y,epsilon = .15,mu = 0,tau2 = 1){

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

  tval = qt(1-epsilon/2,df = N-1)

  ci_help = tval * sqrt(varj + sig2)

  out = list("bounds" = c(thetaj - ci_help,thetaj + ci_help),
             "coverage" = (1-epsilon)*100,
             "data" = Y, "class" = "continuous")

  class(out) = 'pred'

  # return pred object
  return(out)
}
#' Obtain a distance-to-average conformal prediction interval
#'
#' This function computes a conformal prediction region under the distance-from-average
#' non-conformity measure. That is, |a + bz*| <= |ci + di z^*| where i indexes training data.
#'
#' @param Y Observed data vector
#' @param epsilon Prediction error rate
#' @return  pred object
#' @export
dtaPrediction = function(Y,epsilon = .15){

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
  k = floor(epsilon*(N+1))
  int = c(S[k],S[2*N-k+1])

  # if ( k == 0 ) {
  #   stop("ALPHA is set too low.")
  # }

  out = list("bounds" = int, "coverage" = (1-k/(N+1))*100,
             "data" = Y, "class" = "continuous")

  class(out) = 'pred'

  # return pred object
  return(out)
}
#' Obtain a pivot prediction interval
#'
#' This function computes a prediction interval under assumed normality.
#'
#' @param Y Observed data vector
#' @param epsilon Prediction error rate
#' @return  pred object
#' @export
normalPrediction = function(Y,epsilon = .15){

  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    print("Y converted to vector!")
  }

  # standard prediction for single group
  N = length(Y)
  ybar = mean(Y)

  tval = qt(1-epsilon/2,N-1)

  sbar = sd(Y) # use county specific sd

  ci_help = tval * sbar * sqrt(1/N+1)

  out = list("bounds" = c(ybar - ci_help,ybar + ci_help),
             "coverage" = (1-epsilon)*100,
             "data" = Y, "class" = "continuous")

  class(out) = 'pred'

  # return pred object
  return(out)
}


