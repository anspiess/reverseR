lmInfl <- function(model, alpha = 0.05, verbose = TRUE, ...) {
  
  if (class(model) != "lm") stop("Not an 'lm' model!")
  
  ## original parameters
  DATA <- eval(model$model)
  X <- DATA[, 2]
  Y <- DATA[, 1]
  PVAL.old <- summary(model)$coefficients[2, 4]
  SLOPE.old <- coef(model)[2]
  INTER.old <- coef(model)[1]
  SE.old <- summary(model)$coefficients[2, 2]
  wts <- model$weights
  
  ## print standard parameters
  if (verbose) {
    if (PVAL.old <= alpha) cat("\nModel is significant at p = ", PVAL.old, ".\n", sep = "")
    else cat("\nModel is insignificant at p = ", PVAL.old, ".\n", sep = "")
  }
  
  ## define change direction
  if (PVAL.old > alpha) change <- "down" else change <- "up"
  
  ## classical influence measures
  INFL <- influence.measures(model)[[1]]
  
  ## create Leave-One-Out models
  modList <- vector("list", length = length(X))
    
  for (i in 1:length(X)) {
    modList[[i]] <- lm(Y[-i] ~ X[-i], weights = wts[-i], ...)
  }
  
  ## get coefficients, p-values, slopes and s.e.(slope)
  Coef <- sapply(modList, function(x) coef(x))
  PVAL.loo <- sapply(modList, function(x) summary(x)$coefficients[2, 4])
  deltaP <- abs(PVAL.old - PVAL.loo)
  Slope <- sapply(modList, function(x)  summary(x)$coefficients[2, 1])
  Inter <- sapply(modList, function(x)  summary(x)$coefficients[1, 1])
  SE <- sapply(modList, function(x)  summary(x)$coefficients[2, 2])
  deltaSlope <- abs(Slope - SLOPE.old)
  deltaInter <- abs(Inter - INTER.old)
  deltaSE <- abs(SE - SE.old)
  RSQ.loo <- sapply(modList, function(x) summary(x)$r.squared)
  SR <- rstudent(model)
  
  ## define Hadi's measure
  hadi <- function(model) {
    h <- hatvalues(model)
    di <- residuals(model)/sqrt(sum(residuals(model)^2))
    k <- length(coef(model))
    h/(1 - h) + k/(1 - h) * di^2/(1 - di^2) 
  }
  HADI <- hadi(model)
  
  ## return extended influence matrix 
  colnames(INFL)[1] <- "dfb.Inter"
  colnames(INFL)[2] <- "dfb.Slope"
  inflMat <- cbind(Idx = 1:nrow(DATA), DATA[, c(2, 1)], Slope = Coef[2, ], Inter = Coef[1, ], 
                   dSlope = deltaSlope, dInter = deltaInter, INFL, hadi = HADI, sR = SR, looP = PVAL.loo, 
                   dP = deltaP, SE = SE, dSE = deltaSE, Rsq = RSQ.loo)
  rownames(inflMat) <- 1:length(X)
  inflMat <- round(inflMat, 3)
  rawMat <- inflMat
  N <- nrow(inflMat)
  
  ## mark influence measures with *...* based on thresholds
  ## Crit 01: P-value reversal
  P <- inflMat[, 16]
  if (PVAL.old <= 0.05) sel <- which(P > 0.05) else sel <- which(P <= 0.05)
  if (length(sel) > 0) inflMat[sel, 16] <- paste("*", inflMat[sel, 16], "*", sep = "")
  
  ## Crit 02: dfbeta slope
  DFB <- inflMat[, 9]
  sel <- which(abs(DFB) > 2/sqrt(N))
  if (length(sel) > 0) inflMat[sel, 9] <- paste("*", inflMat[sel, 9], "*", sep = "")
  
  ## Crit 03: dffit
  DFF <- inflMat[, 10]
  sel <- which(abs(DFF) > 2 * sqrt(2/N))
  if (length(sel) > 0) inflMat[sel, 10] <- paste("*", inflMat[sel, 10], "*", sep = "")
  
  ## Crit 04: covr
  COVR <- inflMat[, 11]
  sel <- which(abs(COVR - 1) > 3 * 2/N)
  if (length(sel) > 0) inflMat[sel, 11] <- paste("*", inflMat[sel, 11], "*", sep = "")
 
  ## Crit 05: Cook's D
  COOK <- inflMat[, 12]
  sel <- which(COOK > qf(0.5, 2, N - 2))
  if (length(sel) > 0) inflMat[sel, 12] <- paste("*", inflMat[sel, 12], "*", sep = "")
  
  ## Crit 05: Leverage
  HAT <- inflMat[, 13]
  sel <- which(HAT > 2 * 2/N)
  if (length(sel) > 0) inflMat[sel, 13] <- paste("*", inflMat[sel, 13], "*", sep = "")
  
  ## Crit 06: Studentized residuals
  SR <- inflMat[, 15]
  sel <- which(abs(SR) > qt(0.975, N - 2 - 1))
  if (length(sel) > 0) inflMat[sel, 15] <- paste("*", inflMat[sel, 15], "*", sep = "")
 
  ## Crit 06: Hadi's measure
  HADI <- inflMat[, 14]
  sel <- which(HADI > median(HADI, na.rm = TRUE) + 2 * mad(HADI, na.rm = TRUE))
  if (length(sel) > 0) inflMat[sel, 14] <- paste("*", inflMat[sel, 14], "*", sep = "")
  
  ## select all influencers 
  if (change == "up") SEL <- which(PVAL.loo > alpha)
  if (change == "down") SEL <- which(PVAL.loo <= alpha)
  
  if (length(SEL) > 0) {
    if (verbose) cat("Found the following significance reversers, in order of strength:\n")
    selP <- PVAL.loo[SEL]
    
    ## order influencers by decreasing strength
    ordP <- if (change == "up") order(selP, decreasing = TRUE) 
            else order(selP, decreasing = FALSE)
    ordSEL <- SEL[ordP]
    
    ## print ordered influencers
    if (verbose) {
      for (i in ordSEL) cat("  Point #", i, " => p = ", rawMat[i, "looP"], ".\n", sep = "")
    }  
  }  
  else {
    selMat <- NULL
    if (verbose) cat("No significance reversers found!\n")
    ordSEL <- NULL
  }
  
  OUT <- list(origModel = model, finalModels = modList, infl = inflMat, raw = rawMat, sel = ordSEL, alpha = alpha, origP = PVAL.old)
  class(OUT) <- "influencers"
  return(OUT)
}

lmPlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## plot full data
  DATA <- eval(infl$origModel$model)
  X <- DATA[, 2]
  Y <- DATA[, 1]
  plot(X, Y, pch = 16, xlab = colnames(DATA)[2], ylab = colnames(DATA)[1], ...)
  grid()
  
  ## add points
  points(X[infl$sel], Y[infl$sel], pch = 16, cex = 1.5, col = "red3", ...)
  
  ## LOO-models regression line
  if (length(infl$sel) > 0) {
    for (i in infl$sel) {
      abline(infl$finalModels[[i]], lwd = 1, col = "red3", ...)
    }
  }
  
  ## original model regression line
  abline(infl$origModel, lwd = 3, col = "black", ...)
}

pvalPlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## extract leave-one-out p-values
  looP <- infl$raw[, "looP"]
  
  ## plot p-values
  par(mar = c(5, 5, 4, 2))
  plot(1:length(looP), looP, pch = 16, type = "o", las = 0, 
       xlab = "Index", ylab = "P-value", 
       log = "y", ylim = c(0.75 * min(looP, na.rm = TRUE), 1.25 * max(looP, na.rm = TRUE)), 
       las = 1, ...)
  
  points(infl$sel, looP[infl$sel], pch = 16, col = "red3", cex = 1.5)
  grid()
  
  ## add full model p-value and alpha border
  abline(h = infl$origP, col = "blue", lwd = 2, ...)
  abline(h = infl$alpha, col = "red3", lwd = 2, ...)
  legend("top", c("Original model P-value", "Selected alpha"), 
         fill = c("blue", "red3"), horiz = TRUE, inset =-0.12, xpd = TRUE, ...)
}

inflPlot <- function(infl, ...) {
  
  Model <- infl$origModel
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
 
  ## extract influence values
  df_a <- abs(infl$raw[, "dfb.Slope"])
  df_fit <- abs(infl$raw[, "dffit"])
  cov_r <- abs(infl$raw[, "cov.r"])
  cook_D <- infl$raw[, "cook.d"]
  lev <-  infl$raw[, "hat"]
  SR <- abs(infl$raw[, "sR"])
  dP <- infl$raw[, "dP"]
  hadi <- infl$raw[, "hadi"]
  
  ## define sample size
  FITTED <- fitted(Model)
  N <- length(FITTED)
  nPar <- length(coef(Model))
  dF <- df.residual(Model)
  
  ## plot influence values vs. delta p-value
  par(mfrow = c(2, 4))
  par(mar = (c(5, 4, 1, 1)))
  
  ## dfslope
  plot(df_a, dP, pch = 16, xlab = "dfbeta slope", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(df_a[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = 2/sqrt(N), col = "darkgreen", lwd = 2, ...) # 2/sqrt(n)
  grid()
  
  ## dffits
  plot(df_fit, dP, pch = 16, xlab = "dffits", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(df_fit[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = 2 * sqrt(nPar/N), col = "darkgreen", lwd = 2, ...)  # 2 * sqrt(k/n)
  grid()
  
  ## covratio
  plot(abs(cov_r - 1), dP, pch = 16, xlab = "Covratio", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(abs(cov_r - 1)[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = 3 * nPar/N, col = "darkgreen", lwd = 2, ...)  # 1 - 3 * k/n
  grid()
  
  ## Cook's D
  plot(cook_D, dP, pch = 16, xlab = "CookD", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(cook_D[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = qf(0.5, nPar, dF), col = "darkgreen", lwd = 2, ...)  
  grid()
  
  ## leverage
  plot(lev, dP, pch = 16, xlab = "Leverage", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(lev[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = 2 * nPar / N , col = "darkgreen", lwd = 2, ...) # 2 * (k/n)
  grid()
  
  ## studentized residual
  plot(SR, dP, pch = 16, xlab = "Studentized residual", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(SR[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = qt(0.975, dF - 1), col = "darkgreen", lwd = 2, ...) # qt(1-alpha/2, df - 1)
  grid()
  
  ## Hadi's measure
  plot(hadi, dP, pch = 16, xlab = "Hadi's measure", ylab = "delta P-value", 
       las = 1, ylim = c(0, max(dP, na.rm = TRUE)), ...)
  points(hadi[infl$sel], dP[infl$sel], pch = 16, col = "red3", cex = 1.5, ...)
  abline(v = median(hadi, na.rm = TRUE) + 2 * mad(hadi, na.rm = TRUE), col = "darkgreen", lwd = 2, ...) # 
  grid()
  
  ## recover default values
  par(mfrow = c(1, 1))
  par(mar = (c(5.1, 4.1, 4.1, 2.1)))
}

slsePlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## slope and SE of original model
  origSlope <- summary(infl$origModel)$coefficients[2, 1]
  origSE <- summary(infl$origModel)$coefficients[2, 2]
  
  ## extract LOO-slopes and LOO-SEs
  Slope <- infl$infl[, "Slope"]
  SE <- infl$infl[, "SE"]
  
  ## plot both 
  plot(Slope, SE, type = "p", pch = 16, xlab = "Slope",
       ylab = "SE", las = 1)
  points(Slope[infl$sel], SE[infl$sel], type = "p", 
         cex = 1.5, pch = 16, col = "red3")
  grid()
  
  ## insert t-border for p = alpha
  MIN <- if(min(Slope) > 0) 0.5 * min(Slope) else 2 * min(Slope)
  MAX <- if(max(Slope) > 0) 2 * max(Slope) else 0.5 * max(Slope)
  Slope <- seq(MIN, MAX, length.out = 200)
  SE <- qt(1 - 0.5 * infl$alpha, df.residual(infl$origModel))
  lines(Slope, abs(Slope/SE), type = "l", col = "black", lwd = 2)
 
  ## insert original model Slope and SE
  abline(v = origSlope, lwd = 2, col = "blue", lty = 2)
  abline(h = origSE, lwd = 2, col = "darkgreen", lty = 2)
  legend("top", c("Original slope", "Original s.e.", "t-border"), 
         col = c("blue", "darkgreen", "black"), lty = c(3, 3, 1), lwd = 3, 
         xpd = TRUE, inset = -0.25, horiz = TRUE, ...)
}