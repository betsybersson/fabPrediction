#' Plot a `pred` object constructed for a categorical response
#'
#' @param x pred object- a list classified as pred containing objects data and bound
#' @param ... additional parameters passed to the default plot method
#' 
#' @return capability to plot pred object. More details: the command `plot(obj)` 
#' plots the empirical densities of each category. Mass denoted in red indicates
#' inclusion in the prediction set
#'
#' @rdname plot
#' @family pred plots
#' @export
plot.pred = function(x, ...){
  
  if (!inherits(x, "pred")){
    stop("Use only with \"pred\" objects!")
  }
  
  if (x$class=="continuous"){
    plot.range = range(x$data,x$bound)
    plot.range = plot.range + diff(plot.range)/20*c(-1,1)
    
    graphics::stripchart(plot.range,
                         type="l",
                         ...)
    graphics::stripchart(x$data,pch="|",add=T)
    graphics::stripchart(x$bounds,type="l",col="red",lwd=3,add=T)
    text(x$bounds,c(1.05,1.05),labels=round(x$bounds,2),col="red")
  } else if (x$class=="categorical"){
  
    cat.names = names(x$test_stats)
    cols = rep("black",length(x$test_stats))
    cols[which(cat.names%in%x$set)] = "red"
    
    mle = as.table(x$data/sum(x$data))
    names(mle) = cat.names
    
    graphics::plot(mle,
                   ylab = "Empirical Probability Mass",
                   col = cols, 
                   ...)
    ### add code to plot only categories included in set- for space
  }
  invisible()
}