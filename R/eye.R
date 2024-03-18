#' Create Identity Matrix
#'
#' This function returns an NxN identity matrix.
#'
#' @param N dimension of square matrix
#' @return NxN identity matrix
#' @export
eye = function(N){
  return(diag(rep(1,N)))
}