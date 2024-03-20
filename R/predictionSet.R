#' Wrapper to obtain a prediction set for categorical data
#'
#' This function computes a prediction set from a number of methods.
#'
#' @param Y Observed data vector
#' @param method Choice of prediction method. Options include FAB, direct, Bayes.
#' @param alpha Prediction mis-coverage rate
#' @param gamma Dirichlet prior concentration for FAB/Bayes methods
#' @param category_names Category names (optional)
#' 
#' @return pred object containing prediction set and interval coverage
#' 
#' @examples
#' 
#' # obtain example categorical data
#' set.seed(1)
#' prob = rdirichlet(50:1)
#' y = rmultinom(1,15,prob)
#' 
#' fab.set = predictionSet(y,
#'   method = "FAB",
#'   gamma = c(50:1))
#' plot(fab.set)
#' 
#' @export
predictionSet = function(Y,method = "FAB",
                         alpha = .15,
                         gamma = rep(1,length(Y)),
                         category_names = 1:length(Y)){
  
  if (method == "FAB"){
    out = fabCategoricalPrediction(Y,alpha = alpha,gamma = gamma,
                                   category_names = category_names)
  } else if (method == "direct"){
    out = fabCategoricalPrediction(Y,alpha = alpha,gamma = rep(0,length(Y)),
                                   category_names = category_names)
  } else if (method == "Bayes"){
    out = bayesMultinomialPrediction(Y,alpha = alpha,gamma = gamma,
                                     category_names = category_names)
  } else {
    stop(paste0("Error! Method ",method," is not a valid option!"))
  }
  
  # return pred object
  return(out)
}