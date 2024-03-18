#' Obtain a Bayesian prediction interval for categorical data
#'
#' This function computes the Bayesian prediction set for a multinomial conjugate family.
#'
#' @param Y Observed data vector of length K containing counts of observations from each of the K categories
#' @param gamma Dirichlet prior concentration for the K categories
#' @param alpha Prediction mis-coverage rate
#' @param category_names Category names (optional)
#' @return pred object
#' @export
bayesMultinomialPrediction = function(Y,alpha = .15,
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


  K = length(Y)
  N = sum(Y)

  ## test
  pz.out  = array(NA,dim=K)
  names(pz.out) = category_names
  for ( i in 1:K ){

    pz.out[i] = ((Y[i] + gamma[i]) >= (Y + gamma)) %*% (Y + gamma) / (N+sum(gamma))

  }

  alpha.set = names(which(pz.out>alpha))

  out = list("set" = alpha.set, "coverage" = "","test_stats" = pz.out,
             "data" = Y,"class" = "categorical")

  class(out) = 'pred'

  # return pred object
  return(out)
}