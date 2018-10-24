lmMult <- function(model, max = 5, n = 10000, alpha = 0.05, method = c("pearson", "spearman"), verbose = TRUE) {
  ## define correlation function for fast p-value calculation
  ## taken from the 'corr.test' function of the 'psych' package
  corr.test <- function (x, y) 
  {
    if (method == "pearson") r <- cor(x, y, use = "complete.obs", method = "pearson")
    else r <- cor(x, y, use = "complete.obs", method = "spearman")
    n <- t(!is.na(x)) %*% (!is.na(y)); n <- min(n)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- 2 * (1 - pt(abs(t), (n - 2)))
    se <- sqrt((1 - r * r)/(n - 2))
    return(list(p = p, t = t, se = se, r = r))
  }
  
  method <- match.arg(method)
  
  ## extract model data
  DATA <- model.frame(model)
  NR <- nrow(DATA)
  
  ## get original p-value
  origP <- summary(model)$coefficients[2, 4]
  
  ## preallocate 
  sampleMat <- matrix(NA, nrow = max * n, ncol = NR)
  pList <- gList <- vector("list", length = max)
  nRev <- rep(NA, max)
  iter <- 1
  
  ## create combination list
  for (i in 1:max) {
    if (verbose) cat("Sampling ", i, " out of ", NR, " => ", sep = "") 
    
    ## preallocate p-value vector and group vector
    pVec <- rep(NA, n)
    gVec <- rep(i, n)
    
    ## for all n, sample n of N, delete samples and calculate p
    for (j in 1:n) {
     SEL <- sample(1:NR, i, replace = FALSE)
     sampleMat[iter, SEL] <- 1
     X <- DATA[-SEL, 2]
     Y <- DATA[-SEL, 1]
     pVec[j] <- corr.test(X, Y)$p
     iter <- iter + 1
    }
    
    ### calculate number of significance reversers
    if (origP <= alpha) nRev[i] <- round(sum(pVec > alpha)/length(pVec) * 100, 2)
    else nRev[i] <- round(sum(pVec <= alpha)/length(pVec) * 100, 2)
    if (verbose) cat(nRev[i], "% significance reversers.\n", sep = "")
    
    ## store p-values and groups in list
    pList[[i]] <- pVec
    gList[[i]] <- gVec
  }
  
  ## concatenate sample matrix and p-values
  pAllVec <- unlist(pList)
  gAllVec <- unlist(gList)
  sampleMat <- cbind(sampleMat, pval = pAllVec, group = gAllVec)
  sampleMat <- as.data.frame(sampleMat)
  
  ## remove redundant samples
  sampleMat <- unique(sampleMat)
  if (verbose) cat(max * n, "samples gave", nrow(sampleMat), "unique combinations.\n")
  
  names(nRev) <- as.character(1:max)
  class(sampleMat) <- c("data.frame", "multinfl")
  OUT <- list(sample = sampleMat, stat = nRev)
  attr(OUT, "stats") <- c(origP, alpha)
  class(OUT) <- "lmMult"
  return(OUT)
}

multPlot <- function(mult, log = FALSE,  ...) {
  x <- mult
  if (class(x) != "lmMult") stop("object is not a result of 'lmMult' !")
  
  ## get original model p-value and alpha border
  origP <- attr(x, "stats")[1]; alpha <- attr(x, "stats")[2]
  X <- x$sample[, "group"]; Y <- x$sample[, "pval"]
  if (log) {
    Y <- log10(Y)
    alpha <- log10(alpha)
  }
  
  ## define point colors: red => reversers
  COL <- if (origP > alpha) ifelse(Y <= alpha, "#CD000055", "#00000055")
         else ifelse(Y > alpha, "#CD000055", "#00000055")
  
  ## plot in 'stripchart' style
  if (log) YLAB <- "log10(P-value)" else YLAB <- "P-value"
  if (log) YLIM <- NULL else YLIM <- c(0, 1)
  X <- X + rnorm(length(X), 0, 0.05)
  plot(X, Y, pch = 16, col = COL, las = 1, ylim = YLIM,
       xlab = "Left-out samples", ylab = YLAB, cex.lab = 1.5, ...)
  abline(h = alpha, col = "black", lwd = 2, lty = 2, ...)
  if (!log) YPOS <- 1.05 else {
    aT <- axTicks(2)
    Y2 <- par()$usr[4]
    Y1 <- aT[length(aT) - 1]
    YPOS <- Y2 + 0.02 * (Y2 - Y1)
  }
  text(1:length(x$stat), YPOS, paste(x$stat, "%"), pos = 3, offset = 1.5, 
       xpd = TRUE, col = "blue3", srt = 90, ...)
}

