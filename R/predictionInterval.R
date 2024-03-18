#' Wrapper to obtain a prediction interval for continuous data
#'
#' This function computes a prediction interval from a number of methods.
#'
#' @param Y Observed data vector
#' @param method Choice of prediction method. Options include FAB, DTA, direct, Bayes.
#' @param alpha Prediction error rate
#' @param mu Prior expected mean of the population mean
#' @param tau2 Prior expected variance of the population mean
#' @return pred object containing prediction interval bounds and interval coverage
#' @export
predictionInterval = function(Y,method = "FAB",
                              alpha = .15,
                              mu = 0,tau2 = 1){
  
  if (method == "FAB"){
    out = fabContinuousPrediction(Y,alpha = alpha,mu = mu,tau2 = tau2)
  } else if (method == "DTA"){
    out = dtaPrediction(Y,alpha = alpha)
  } else if (method == "direct"){
    out = normalPrediction(Y,alpha = alpha)
  } else if (method == "Bayes"){
    out = bayesNormalPrediction(Y,alpha = alpha,mu = mu,tau2 = tau2)
  } else {
    stop(paste0("Error! Method ",method," is not a valid option!"))
  }
  
  # return pred object
  return(out)
}