#' Obtain empirical Bayesian estimates for conjugate normal spatial Fay-Herriot model
#'
#' This function returns plug-in values for a conjugate
#' normal spatial Fay-Herriot model.
#'
#' @param Y Data vector
#' @param group Group membership of each entry in Y
#' @param W Adjacency matrix
#' @param X Group-level covariates
#' @return plug-in values of spatial Fay-Herriot model
#' @export
plugin_values = function(Y,group,W = NA,X = NA){

  if (!is.factor(group)){
    group = as.factor(group)
  }

  ybar = tapply(Y,group,mean)
  N = table(group)
  ss = (N-1)*c(tapply(Y,group,var))

  ## obtain mles of s2
  mll = function(XX){
    a = XX[1]
    b = XX[2]

    ap = (a+N-1)/2
    bp = (b+ss)/2

    out = -sum( ( (a/2)*log(b/2)-lgamma(a/2) ) - ( ap*log(bp)-lgamma(ap) ) )

    return(out)
  }
  init = c(1,1)
  ab_out = optim(init, mll)$par

  ## posterior modes for each s2
  a = ab_out[1]
  b = ab_out[2]

  ap = (a+N-1)/2
  bp = (b+ss)/2

  # est of var
  s2 = bp/(ap+1)

  if ( any(is.na(W))|all(W==eye(ncol(W))) ){ # if W is empty or the identity- don't use spatial method

    q = ncol(X)

    ## obtain mle of mu, t2 with plugin s2
    mll2 = function(XX){
      tau2 = XX[1]
      mu = XX[-1]

      out = -sum(dnorm(ybar,X%*%mu,sqrt(s2*(1/N+exp(tau2))),log=TRUE) )

      return(out)
    }
    init = c(log(var(ybar)), rep(mean(ybar),q) )
    mutau_out = optim(init,mll2)$par
    mutau_out[1] = exp(mutau_out[1])

    out = list(mu=mutau_out[-1],
               tau2=mutau_out[1])

  } else {

    ## spatial ests of mu, t2, rho, theta using plugin s2

    ## organize data frame
    df = data.frame(cbind(c(ybar),c(s2/N),X))
    colnames(df)[1:2] = c("Y","Var")

    resultML = sae::eblupSFH(Y ~ . - Var - 1, Var, as.data.frame(W),
                        method="ML",data = df)

    theta.eb = resultML$eblup
    mu.eb = resultML$fit$estcoef[1]
    rho.eb = resultML$fit$spatialcorr
    tau2.eb = resultML$fit$refvar

    out = list(mu = mu.eb,
               tau2 = tau2.eb,
               rho = rho.eb,
               theta = theta.eb,
               s20 = b/(a+1))
  }

  return(out)

}

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
EB_values = function(j,Y,group,W = NA,X = NA){

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

    out = plugin_values(Y[j!=group],rep(1:(J-1),ns[-j]),
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

    out = plugin_values(Y[j!=group],rep(1:(J-1),ns[-j]),
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

##################################



