#' Obtain Hessian of the marginal Dirichlet-multinomial likelihood
#'
#'
#' @param D matrix (JxK) of counts; each row is a sample from a MN distribution with K categories
#' @param gamma current value of prior concentration parameter
#' @param Nj sample sizes of the J groups
#' @param K number of categories
#' @return Hessian
#' @export
polyaHessian= function(D, gamma,
                       Nj = rowSums(D),
                       K = ncol(D)){
  ## helpers
  gamma0 = sum(gamma)
  
  z = sum(trigamma(gamma0) - 
            trigamma(gamma0 + Nj))
  
  q.diag = unlist(parallel::mclapply(1:K,function(k)
    sum(trigamma(D[,k] + gamma[k]) -
          trigamma(gamma[k]))
  ))
  Q = diag(q.diag)
  
  H = Q + z
  
  return(list("q.diag" = q.diag,"z" = z,"H" = H))
}