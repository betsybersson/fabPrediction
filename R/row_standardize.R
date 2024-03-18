#' Row standardize a matrix
#'
#'
#' @param W matrix
#' @return row-standardized matrix
#' @export
row_standardize = function(W){
  return(diag(1/rowSums(W)) %*% W)
}