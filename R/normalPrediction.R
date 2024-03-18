#' Obtain a pivot prediction interval
#'
#' This function computes a prediction interval under assumed normality.
#'
#' @param Y Observed data vector
#' @param alpha Prediction error rate
#' @return  pred object
#' @export
normalPrediction = function(Y,alpha = .15){

  if (!is.vector(Y)){
    Y = unlist(as.vector(Y))
    Y = unname(Y)
    print("Y converted to vector!")
  }

  # standard prediction for single group
  N = length(Y)
  ybar = mean(Y)

  tval = qt(1-alpha/2,N-1)

  sbar = sd(Y) # use county specific sd

  ci_help = tval * sbar * sqrt(1/N+1)

  out = list("bounds" = c(ybar - ci_help,ybar + ci_help),
             "coverage" = (1-alpha)*100,
             "data" = Y, "class" = "continuous")

  class(out) = 'pred'

  # return pred object
  return(out)
}