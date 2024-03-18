#' Obtain gradient of the marginal Dirichlet-multinomial likelihood
#'
#'
#' @param D matrix (JxK) of counts; each row is a sample from a MN distribution with K categories
#' @param gamma current value of prior concentration parameter
#' @param Nj sample sizes of the J groups
#' @param K number of categories
#' @return gradient
#' @export
polyaGradient = function(D, gamma,
                         Nj = rowSums(D),
                         K = ncol(D)){
  ## helpers
  gamma0 = sum(gamma)
  
  g_k = unlist(
    parallel::mclapply(1:K,function(k) sum(digamma(gamma0) -
                                             digamma(gamma0 + Nj) +
                                             digamma(D[,k] + gamma[k]) -
                                             digamma(gamma[k])))
  )
  
  return(g_k)
}