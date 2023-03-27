#' Create Identity Matrix
#'
#' This function returns an NxN identity matrix.
#'
#' @param N dimension of square matrix
#' @return NxN identity matrix
#' @export
eye = function(N){
  diag(rep(1,N))
}

#' Plot a `pred` object
#'
#'
#' @param obj pred object- a list classified as pred containing objects data and bound
#' @param main plot title
#' @param xlab xaxis title
#' @return capability to plot pred object. More details: the command `plot(obj)` returns a line plot with tick marks denoting observed data points and a red band denoting the prediction interval
#' @export
plot.pred = function(obj,
                     main="",
                     xlab=""){

if (obj$class=="continuous"){
  plot.range = range(obj$data,obj$bound)
  plot.range = plot.range + diff(plot.range)/20*c(-1,1)

  stripchart(plot.range,type="l",main=main,xlab=xlab)
  stripchart(obj$data,pch="|",add=T)
  stripchart(obj$bounds,type="l",col="red",lwd=3,add=T)
  text(obj$bounds,c(1.05,1.05),labels=round(obj$bounds,2),col="red")
}

}
