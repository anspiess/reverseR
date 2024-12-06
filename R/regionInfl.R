regionInfl <- function(model, div.x = 20, div.y = 20, grid = TRUE, pred.int = TRUE,
  crit = c("P", "dfb.Slope", "dffit", "cov.r", "cook.d", "hat", "hadi", "sR", "cdr", "Si"),
  cex.grid = 0.5, alpha = 0.05, xlim = NULL, ylim = NULL, ...)
{
  if (is.null(xlim) | is.null(ylim)) stop("'xlim' and 'ylim' must be defined!")
  crit <- match.arg(crit)
  DAT <- model.frame(model)
  X <- DAT[, 2]; Y <- DAT[, 1]
  LM <- lm(Y ~ X)
  
  plot(X, Y, type = "n", xlim = xlim, ylim = ylim, mgp = c(2, 0.5, 0), ...)
  USR <- par()$usr

  seqX <- seq(xlim[1], xlim[2], length.out = div.x)
  seqY <- rev(seq(ylim[1], ylim[2], length.out = div.y))
  
  if (grid) {
    abline(v = seqX, col = "grey80")
    abline(h = seqY, col = "grey80") 
  }
  
  Grid <- expand.grid(seqX, seqY); Grid <- cbind(Grid, "crit" = NA_real_)
  colnames(Grid) <- c("x", "y", "crit")
  pb <- txtProgressBar(min = 1, max = nrow(Grid), initial = 0, style = 3, char = "-")
  
  for (i in 1:nrow(Grid)) {
    setTxtProgressBar(pb, i)
    addX <- Grid[i, 1]; addY <- Grid[i, 2]
    X2 <- c(X, addX); Y2 <- c(Y, addY) 
    LEN <- length(X2)
    addModel <- lm(Y2 ~ X2)
    
    if (crit == "P") {
      PVAL <- summary(addModel)$coefficients[2, 4]
      if (PVAL <= alpha) COL <- "#006400" else COL <- "#FF8C00"
      if (PVAL <= alpha) Grid[i, 3] <- 1 else Grid[i, 3] <- 0
    } else {
      addInfl <- lmInfl(addModel, verbose = FALSE, ...)
      if (grepl("\\*", addInfl$infl[LEN, crit])) COL <- "#FF8C00" else COL <- "#006400"
      if (grepl("\\*", addInfl$infl[LEN, crit])) Grid[i, 3] <- 1 else Grid[i, 3] <- 0
    }
    
    points(addX, addY, pch = 16, cex = cex.grid, col = COL)
  }
  
  points(X, Y, pch = 16, cex = 2)
  abline(LM, col = "black", lwd = 2, lty = 1)
  
  PRED <- predict(LM, interval = "prediction", newdata = data.frame(X = seqX), level = 1 - alpha)
  if (pred.int) {
    lines(seqX, PRED[, 2], lty = 2, lwd = 2, col = "black")
    lines(seqX, PRED[, 3], lty = 2, lwd = 2, col = "black")
  }
  
  return(grid = Grid)
}
