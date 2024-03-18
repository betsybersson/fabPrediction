#' Obtain inital guess of MLE of the marginal Dirichlet-multinomial likelihood
#'
#' Method of moment matching to obtain an initial guess of the MLE, as in Minka (2000).
#'
#' @param D matrix (JxK) of counts; each row is a sample from a MN distribution with K categories
#' @return Hessian
#' @export
initMoM = function(D){
  p = t(apply(D,1,function(l)l/sum(l)))
  exp.p = colMeans(p)
  exp2.p = colMeans(p*p)
  ok = exp.p>0
  sum.as = (exp.p[ok] - exp2.p[ok]) / (exp2.p[ok] - exp.p[ok]^2)
  init = exp.p * median(sum.as)
  return(init)
}