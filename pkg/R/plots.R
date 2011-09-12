#' Display a heatmap of factor loadings (note that this function requires ggplot2 package)
#' @aliases plot_loadings
#' @param x an sbfac model object or a (loadings) matrix
#' @param xlabel x-axis label
#' @param ylabel y-axis label
#' @param color character vector of colors (default Bl-Wh-Rd)
#' @param sorting a permutaion of 1:P (where P is the number of variables) providing sort order
#' @param type Either "color" (heatmap) or "line" (B&W alternative)
#' for the rows of the loadings matrix
#' @param scale.name label for the legend
#' @export

plot_loadings <- function(x, xlabel=NA, ylabel=NA, 
                    color=NA, sorting=NA, scale.name="Loading", type="color") {
  if (class(x)=="bfa") {
    loadings = x$post.loadings.mean
    rownames(loadings) = x$varlabel
  } else {
    loadings = x
  }
  if (is.na(color)) {
    color=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
            "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
  }
  
	ldf = melt(loadings)
	colnames(ldf) = c("Group", "Factor", "value")
  if (!any(is.na(sorting))) ldf$Group = factor(ldf$Group, levels=rownames(loadings)[sorting])
  ldf$Factor = as.factor(ldf$Factor)
  
  if (type=="color") {
	breaks = seq(-1,1,by = 0.2)
	cl = colorRampPalette(color)(21)
  	lim = c(-1.0, 1)
  	p = ggplot(ldf, aes(x=Factor, y=Group, fill=value))
  	p = p + scale_fill_gradientn(scale.name, colour=cl, limits = lim, breaks=breaks) 
    p = p + geom_tile()
  }
  
  if (type=="line") {
    k = x$K
    numgp = x$P
    
    ldf$numerfac = as.numeric(ldf$Factor)
    ldf$xend     = ldf$numerfac+0.45*ldf$value
    
    xi = c(0.55, (1:k)+0.45, (1:(k-1))+0.55)
    
    p = ggplot(ldf, aes(x=numerfac, y=Group))
    p = p + geom_vline(xintercept=1:k, colour="gray70")
    p = p + geom_hline(aes(yintercept=Group), colour="gray86")
    p = p + geom_vline(xintercept=xi, colour="gray70", linetype=2)
    p = p + geom_segment(aes(xend = xend, yend=Group), size=2)
    p = p + opts(panel.grid.y.major = theme_line(colour = 'gray60', linetype = 1))
    p = p + geom_rect(aes(xmin = (0:k)+0.485, xmax=(0:k)+0.515, ymin=rep(0,k), ymax=rep(numgp+1, k)), 
                      fill="gray80", color="gray80")
    p = p+scale_y_discrete(limits=levels(ldf$Group))+xlim(c(0.5, k+0.45))
    p = p + theme_bw() + scale_x_continuous("Factor", breaks=1:k)
  }
  
  p = p+ ylab("")+xlab("")
	if (!is.na(xlabel)) p = p+xlab(xlabel)
	if (!is.na(ylabel)) p = p+ylab(ylabel) 
  return(p)

}

#' Display a biplot
#' @param x A bfa object
#' @param factors Numeric/integer vector giving indices of the factors to plot
#' @param ... Additional arguments to biplot; see \code{?biplot}
#' @method biplot bfa
#' @return Shows a biplot
#' @export
biplot.bfa <- function(x, factors=c(1,2), ...) {
  call_args = list(...)
  if(is.null(call_args$xlabs)) call_args$xlabs = x$obslabel
  if(is.null(call_args$ylabs)) call_args$ylabs = x$varlabel
  call_args$x = t(x$post.scores.mean)[,factors]
  call_args$y = x$post.loadings.mean[,factors]
  do.call(biplot, call_args)
}