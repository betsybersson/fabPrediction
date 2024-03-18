#' Obtain empirical Bayesian estimates for group j
#'
#' This function returns empirical Bayesian estimates for a specified group from
#' the conjugate normal spatial Fay-Herriot model.
#'
#' @param j Obtain EB values for group in index j- numeric value in group
#' @param Y Data vector
#' @param group index vecter of the same lenght as Y
#' @param W Non-standardized adjacency matrix
#' @param X Group-level covariates
#' @return empirical Bayesian estimates of population mean and it's variance
#' 
#' @export
fayHerriotEB = function(j,Y,group,W = NA,X = NA){

  # dimensions
  ns = table(group)
  J = length(ns)

  ## if no covariates entered, use covariates
  if (any(is.na(X))){
    if(!all(is.na(X))){
      warning("NA in covariate input- just using intercept")
    }
    X = matrix(1,ncol=1,nrow=J)
  }

  ## if X is vector, convert to matrix.
  if (is.vector(X)) {
    X = matrix(X,ncol=1)
  }

  if ( any(is.na(W))|all(W==eye(ncol(W))) ){ # if W is empty or the identity- don't use spatial method

    out = pluginValues(Y[j!=group],rep(1:(J-1),ns[-j]),
                        X[-j,,drop=F])

    # extract params
    mu = unlist(out$mu)
    tau2 = c(out$tau2)

    # output
    MU = c(X[j,,drop=F] %*% mu)
    TAU2 = tau2

    out = list("mu" = MU,"tau2" = TAU2)

  } else {

    if (all(rowSums(W)==1)){
      warning("Please input the non-row-standardized adjacency matrix!")
    }

    W.stand = row_standardize(W) #row standardize

    out = pluginValues(Y[j!=group],rep(1:(J-1),ns[-j]),
                        W.stand[-j,-j],X[-j,,drop=F])

    # extract params
    mu = unlist(out$mu)
    rho = c(out$rho)
    omega2 = c(out$tau2)
    theta = c(out$theta)
    s2 = c(out$s20)

    # helper
    G = omega2 * solve((eye(J)-rho*t(W.stand))%*%(eye(J)-rho*(W.stand)))
    R = G[j,-j] %*% solve(G[-j,-j])

    # output
    MU = c(X[j,,drop=F] %*% mu + R %*% (theta - X[-j,,drop=F] %*% mu))
    TAU2 = c(G[j,j] - R %*% G[-j,j])/s2

    out = list("mu" = MU,"tau2" = TAU2)
  }

  return(out)

}