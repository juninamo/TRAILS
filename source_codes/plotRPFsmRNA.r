
plotRPFsmRNA <- function(TE, sample, alpha, cex, color_TE=FALSE, 
                   removeZero=TRUE, log2=TRUE, breaks.length=50, ...){
  if(!is.list(TE)){
    stop("TE must be output of translationalEfficiency.")
  }
  if(!any(c("RPFs", "mRNA", "TE") %in% names(TE))){
    stop("TE must be output of translationalEfficiency.")
  }
  if(missing(sample)){
    stop("sample is required.")
  }
  norm01 = function(x){(x-min(x))/(max(x)-min(x))}
  rbPal <- colorRampPalette(c('red','blue'))
  
  mRNA <- TE$mRNA
  RPFs <- TE$RPFs
  TE <- TE$TE
  if(!is.numeric(sample)){
    sample <- which(colnames(TE) %in% sample)
  }
  if(length(sample)>1){
    sample <- sample[1]
    message("Only first sample will be plotted.")
  }
  TE <- TE[, sample]
  RPFs <- RPFs[, sample]
  mRNA <- mRNA[, sample]
  if(removeZero){
    keep <- RPFs>0 & mRNA>0
    if(sum(keep)<1){
      stop("No data available for plotting.")
    }
    mRNA <- mRNA[keep]
    RPFs <- RPFs[keep]
    TE <- TE[keep]
  }
  if(log2){
    mRNA <- log2(mRNA)
    RPFs <- log2(RPFs)
    TE <- log2(TE)
  }
  test = cor.test(mRNA,RPFs, na.rm = TRUE)
  opar <- par(fig=c(0, .75, 0, .75), new=FALSE, mar=c(5.1, 4.1, 0, 0))
  on.exit(par(opar))
  dots <- list(...)
  args <- dots
  args$x <- RPFs
  args$y <- mRNA
  if (!color_TE){
    args$col <- alpha
  } else {
    args$col <- rbPal(10)[as.numeric(cut(norm01(TE),breaks = 10))]
  }
  if(length(args$xlab)==0) args$xlab <- "RPFs level"
  if(length(args$ylab)==0) args$ylab <- "mRNA level"
  do.call(plot, args)
  par(fig=c(0, .75, 0, .75), new=TRUE, mar=c(5.1, 4.1, 0, 0))
  args <- dots
  args$x <- quantile(RPFs,probs = seq(0, 1, 0.0001), na.rm = TRUE)["95%"]
  args$y <- quantile(mRNA,probs = seq(0, 1, 0.0001), na.rm = TRUE)["99.999%"]
  args$col <- "black"
  args$cex <- cex
  args$labels <- paste0(ifelse(test$p.value == 0, "p<2.2e-16",format(test$p.value, digit=2)),", R=",format(test$estimate, digit=2))
  do.call(text, args)
  ylim <- par("usr")[3:4]
  xlim <- par("usr")[1:2]
  par(fig=c(.75, 1, 0, .75), new=TRUE, mar=c(5.1, 0, 4.1, 5.1))
  yhist <- hist(mRNA, breaks=seq(ylim[1], ylim[2], length.out = breaks.length),
                plot=FALSE)
  args <- dots
  args$height <- yhist$density
  args$axes <- FALSE
  args$space <- 0
  args$horiz <- TRUE
  args$cex <- NULL
  do.call(barplot, args)
  par(fig=c(0, .75, 0.75, 1), new=TRUE, mar=c(0, 4.1, 4.1, 0))
  xhist <- hist(RPFs, breaks=seq(xlim[1], xlim[2], length.out = breaks.length),
                plot=FALSE)
  args <- dots
  args$height <- xhist$density
  args$axes <- FALSE
  args$space <- 0
  args$horiz <- FALSE
  args$cex <- NULL
  do.call(barplot, args)
}