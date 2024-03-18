#' Generate a random sample from a Dirichlet distribution
#'
#'
#' @param gamma Prior concentration vector of length K
#' @return a vector of length K that is a random sample from a Dirichlet distribution
#' @export
rdirichlet = function(gamma){
  Y = sapply(gamma,function(j)rgamma(1,j,1))
  
  return(Y/sum(Y))
}