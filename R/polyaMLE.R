#' Obtain MLE of marginal Dirichlet-multinomial likelihood
#'
#' This function retuns the MLE of the prior concentration from a marginal Dirichlet-multinomial
#' likelihood. Default method iterates a Newton-Raphson algorithm until convergence.
#'
#' @param D matrix (JxK) of counts; each row is a sample from a MN distribution with K categories
#' @param init If NA, use method moment matching procedure to obtain good init values
#' @param method "Newton_Raphson", "fixed_point", "separate", "precision_only"
#' @param epsilon convergence diagnostic
#' @param print_progress if TRUE, print progress to screen
#' @return mle of prior concentration from marginal Dirichlet-multinomial likelihood
#' @export
polyaMLE = function(D, init = NA, method = "Newton_Raphson",
                    epsilon = .0001,
                    print_progress = FALSE) {
  
  ## if a column is all 0s, remove it from analysis and set prior to 0
  zeros.ind = which(colSums(D)==0)
  if (length(zeros.ind)>0){
    not.zeros.ind = c(1:ncol(D))[-zeros.ind]
    D.orig = D
    D = D[,-zeros.ind]
  }
  
  if (is.na(init)){
    init = initMoM(D) 
  }
  
  # grab helpers
  Nj = rowSums(D)
  J = nrow(D)
  K = ncol(D)
  
  # some initialization steps
  alpha_t = init
  if (method == "separate"| method == "precision_only"){
    theta0_t = colMeans(row_standardize(D))
    alpha0_t = mean(rowSums(D))
    nu.fn = function(nn,alpha.star){
      sum(alpha.star * (digamma(nn + alpha.star) - digamma(alpha.star)))
    }
  }
  converged = FALSE
  # iterate until convergence
  while (!converged) {
    
    if (method == "Newton_Raphson"){
      # get gradient and inverse Hessian
      g = polyaGradient(D,alpha_t,Nj,K)
      H = polyaHessian(D,alpha_t,Nj,K)$H
      H.inv = solve(H) # should be able to reduce comp time a bit more than this
      
      # update alpha 
      alpha_tP1 = alpha_t - c(H.inv %*% g)
    } else if (method == "fixed_point") {
      # get helper
      alpha0t = sum(alpha_t)
      gprime.k = unlist(parallel::mclapply(1:K,function(k) sum(digamma(D[,k]+alpha_t[k]) -
                                                                 digamma(alpha_t[k]))/
                                             sum(digamma(Nj + alpha0t) - digamma(alpha0t))
      ))
      
      # update alpha
      alpha_tP1 = alpha_t*gprime.k
    } else if (method == "separate"){
      # update precision- alpha0
      gprime.k = sum(unlist(parallel::mclapply(1:K,function(k) 
        sum(theta0_t[k] * digamma(D[,k]+ alpha0_t * theta0_t[k]) -
              theta0_t[k] * digamma(alpha0_t * theta0_t[k])))))/
        sum(digamma(Nj + alpha0_t) - digamma(alpha0_t))
      # update alpha0
      alpha0_tP1 = alpha0_t*gprime.k
      
      ## update theta0
      theta0k =  unlist(parallel::mclapply(1:K,function(k)nu.fn(D[,k],alpha0_tP1*theta0_t[k])))
      theta0_tP1 = theta0k/sum(theta0k)
      
      # compute alphat for convergence check and output
      alpha_tP1 = theta0_tP1*alpha0_tP1
      
      # save updated values
      alpha0_t = alpha0_tP1
      theta_t = theta0_tP1
      
    }else if (method == "precision_only"){
      # update precision- alpha0
      gprime.k = sum(unlist(parallel::mclapply(1:K,function(k) 
        sum(theta0_t[k] * digamma(D[,k]+ alpha0_t * theta0_t[k]) -
              theta0_t[k] * digamma(alpha0_t * theta0_t[k])))))/
        sum(digamma(Nj + alpha0_t) - digamma(alpha0_t))
      # update alpha0
      alpha0_tP1 = alpha0_t*gprime.k
      
      # compute alphat for convergence check and output
      alpha_tP1 = theta0_t*alpha0_tP1
      
      # save updated values
      alpha0_t = alpha0_tP1
    }
    
    
    # check convergence
    converged = ifelse(all(abs(alpha_tP1-alpha_t)<epsilon),TRUE,FALSE)
    alpha_t = ifelse(alpha_tP1<0,1e-10,alpha_tP1) # enforce bounded below by 0 constraint
    
    if (print_progress){
      print(paste0(c("alpha:",round(alpha_t,2)), collapse= " "))
    }
    
  }
  
  # if some species are all 0, put prior of 0 back in
  if (length(zeros.ind)>0){
    alpha.out = array(0,dim=ncol(D.orig))
    alpha.out[not.zeros.ind] = alpha_t
    alpha_t = alpha.out
  }
  
  return(alpha_t)
  
}