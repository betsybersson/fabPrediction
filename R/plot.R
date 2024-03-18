#' Plot a `pred` object
#'
#'
#' @param obj pred object- a list classified as pred containing objects data and bound
#' @param plot.all (TO BE ADDED) T displays all categories when plotting a categorical prediction set.
#' @return capability to plot pred object. More details: the command `plot(obj)` returns a line plot with tick marks denoting observed data points and a red band denoting the prediction interval
#' @export
plot.pred = function(obj,
                     main="",
                     xlab="",
                     cex.axis=1){
  
  if (obj$class=="continuous"){
    plot.range = range(obj$data,obj$bound)
    plot.range = plot.range + diff(plot.range)/20*c(-1,1)
  
    stripchart(plot.range,type="l",main=main,xlab=xlab)
    stripchart(obj$data,pch="|",add=T)
    stripchart(obj$bounds,type="l",col="red",lwd=3,add=T)
    text(obj$bounds,c(1.05,1.05),labels=round(obj$bounds,2),col="red")
  } else if (obj$class=="categorical"){
    cat.names = names(obj$test_stats)
    cols = rep("black",length(obj$test_stats))
    cols[which(cat.names%in%obj$set)] = "red"
    
    mle = as.table(obj$data/sum(obj$data))
    names(mle) = cat.names
  
  
    plot.all = T
    if (plot.all == T){
      plot(mle,
           ylab = "Empirical Probability Mass",
           col = cols,
           xlab = xlab,
           main = main,
           cex.axis = cex.axis)
    } else {
      ### add code to plot only categories included in set- for space
    }
    
  }
  invisible()
}