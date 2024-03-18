#' Obtain a FAB conformal prediction interval for categorical data
#'
#' This function computes a FAB conformal prediction set as described in
#' Bersson and Hoff 2023.
#'
#' @param Y Observed data vector of length K containing counts of observations from each of the K categories
#' @param gamma Dirichlet prior concentration for the K categories
#' @param alpha Prediction mis-coverage rate
#' @param category_names Category names (optional)
#' @return pred object
#' @export
fabCategoricalPrediction = function(Y, alpha = .15,
                                    gamma = rep(1,length(Y)),
                                    category_names = 1:length(Y)){
  
  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    print("Y converted to vector!")
  }
  
  K = length(Y)
  if ( K != length(gamma)){
    stop("gamma must be the same length as data count vector Y")
  }
  if ( K != length(category_names)){
    stop("category_names must be the same length as data count vector Y")
  }
  N = sum(Y)
  
  
  # for each category, obtain xl + gammal
  c_i = Y + gamma
  
  # for each category, obtain c_nPLUS1 in that case
  c_nPLUS1 = Y + 1 + gamma
  
  ## test
  pz.out = array(NA,dim=K)
  names(pz.out) = category_names
  for ( k in 1:K ){
    
    pz = c(c_nPLUS1[k] >= c_i) %*% Y + 1
    # normalize
    pz = pz/(N+1)
    
    pz.out[k] = pz
    
  }
  
  alpha.set = names(which(pz.out>alpha))
  
  out = list("set" = alpha.set, "coverage" = "", "test_stats" = pz.out,
             "data" = Y, "class" = "categorical")
  
  class(out) = 'pred'
  
  # return pred object
  return(out)
}