inflPlot <- function(infl, measure = NULL, ...) {
  
  ## check for correct class
  if (as.character(class(infl)) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## check for correct measure
  CN <- colnames(infl$infl)[8:22]
  if (!is.null(measure) && !(measure %in% CN)) stop(cat("'measure' must be any of ", CN, "!\n"))
  
  ## plot full data
  DATA <- eval(infl$origModel$model)
  X <- DATA[, 2]
  Y <- DATA[, 1]
  par(mar = c(5.1, 5.1, 4.1, 1.1))
  plot(X, Y, pch = 16, xlab = colnames(DATA)[2], ylab = colnames(DATA)[1], panel.first = grid(), cex.lab = 2, cex.axis = 1.5, cex = 1.5, col = "gray30", ...)
  abline(infl$origModel, lwd = 3, col = "black", ...)
 
  ## add points
  if (!is.null(infl$sel)) {
    SEL <- infl$sel[1]
    points(X[infl$sel], Y[infl$sel], pch = 16, cex = 2, col = "darkred", ...)
    abline(infl$finalModels[[SEL]], col = "darkred", lwd = 3, ...)
  }
  
  ## legend
  legend("topleft", legend = c("orig. trend", "reverser trend"), lty = 1, lwd = 3, col = c("black", "darkred"), horiz = TRUE, xpd = TRUE, inset = c(0, -0.08))
}

pvalPlot <- function(infl, ...) {
  
  ## check for correct class
  if (as.character(class(infl)) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## extract leave-one-out p-values
  looP <- infl$raw[, "looP"]
  
  ## plot p-values
  par(mar = c(5.1, 5.1, 4.1, 1.1))
  plot(1:length(looP), looP, pch = 16, type = "o", xlab = "Index", ylab = "log(P-value)", col = "black",
       ylim = c(0.75 * min(looP, na.rm = TRUE), 1.25 * max(looP, na.rm = TRUE)), log = "y", cex.lab = 2, cex.axis = 1.5, panel.first = grid(), ...)
  
  ## add full model p-value and alpha border
  abline(h = infl$alpha, col = "darkgreen", lwd = 2, ...)
  abline(h = infl$origP, col = "black", lwd = 2, lty = 2, ...)
  
  ## color reversers
  if (!is.null(infl$sel)) {
    SEL <- infl$sel
    points(SEL, looP[SEL], pch = 16, col = "darkred", cex = 2)
  }
  
  ## legend
  legend("topleft", legend = c("alpha", "orig. P-value"), lty = c(1, 2), lwd = 3, col = c("darkgreen", "black"), horiz = TRUE, xpd = TRUE, inset = c(0, -0.08))
}
