lmInfl <- function(model, alpha = 0.05, cutoff = c("BKW", "R"), verbose = TRUE, ...) {
  
  if (as.character(class(model)) != "lm") stop("Not an 'lm' model!")
  cutoff <- match.arg(cutoff)
  
  ## original parameters
  DATA <- eval(model$model)
  X <- DATA[, 2]
  Y <- DATA[, 1]
  PVAL.old <- summary(model)$coefficients[2, 4]
  TVAL.old <- summary(model)$coefficients[2, 3]
  SLOPE.old <- coef(model)[2]
  INTER.old <- coef(model)[1]
  SE.old <- summary(model)$coefficients[2, 2]
  RSQ.old <- summary(model)$r.squared
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
    y <- Y[-i]; x <- X[-i]
    modList[[i]] <- lm(y ~ x, weights = wts[-i], ...)
  }
  
  ## get coefficients, p-values, slopes and s.e.(slope)
  Coef <- sapply(modList, function(x) coef(x))
  PVAL.loo <- sapply(modList, function(x) summary(x)$coefficients[2, 4])
  Dfstat <- TVAL.old - sapply(modList, function(x) summary(x)$coefficients[2, 3])
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
  
  ## define Coefficient of Determination Ratio (CDR)
  CDR <- RSQ.loo/RSQ.old
  
  ## define Pena's measure
  residMat <- matrix(NA_real_, ncol = nrow(DATA), nrow = nrow(DATA))
  fitted0 <- fitted(model)
  for (i in 1:length(modList)) {
    residMat[i, ] <- fitted0 - predict(modList[[i]], newdata = data.frame(x = X))
  }
  Si.sum <- colSums(residMat^2, na.rm = TRUE)
  Si.var <- sum(residuals(model)^2)/df.residual(model) * hatvalues(model)
  Si <- Si.sum/(2 * Si.var)
  
  ## return extended influence matrix 
  colnames(INFL)[1] <- "dfb.Inter"
  colnames(INFL)[2] <- "dfb.Slope"
  inflMat <- cbind(Idx = 1:nrow(DATA), DATA[, c(2, 1)], Slope = Coef[2, ], Inter = Coef[1, ], 
                   dSlope = deltaSlope, dInter = deltaInter, INFL, hadi = HADI, sR = SR, dfstat = Dfstat, looP = PVAL.loo, 
                   dP = deltaP, SE = SE, dSE = deltaSE, Rsq = RSQ.loo, cdr = CDR, Si = Si)
  rownames(inflMat) <- 1:length(X)
  inflMat <- round(inflMat, 3)
  rawMat <- inflMat
  N <- nrow(inflMat)
  
  ## mark influence measures with *...* based on thresholds
  ## Crit 01: P-value reversal
  P <- inflMat[, 17]
  if (PVAL.old <= 0.05) sel <- which(P > 0.05) else sel <- which(P <= 0.05)
  if (length(sel) > 0) inflMat[sel, 17] <- paste("*", inflMat[sel, 17], "*", sep = "")
  
  ## Crit 02: dfbeta slope
  DFB <- inflMat[, 9]
  if (cutoff == "BKW") sel <- which(abs(DFB) > 2/sqrt(N)) else sel <- which(abs(DFB) > 1)
  if (length(sel) > 0) inflMat[sel, 9] <- paste("*", inflMat[sel, 9], "*", sep = "")
  
  ## Crit 03: dffit
  DFF <- inflMat[, 10]
  if (cutoff == "BKW") sel <- which(abs(DFF) > 2 * sqrt(2/N)) else sel <- which(abs(DFF) > 3 * sqrt(2/(N - 2)))
  if (length(sel) > 0) inflMat[sel, 10] <- paste("*", inflMat[sel, 10], "*", sep = "")
  
  ## Crit 04: covr
  COVR <- inflMat[, 11]
  if (cutoff == "BKW") sel <- which(abs(COVR - 1) > 3 * 2/N) else sel <- which(abs(COVR - 1) > 3 * 2/(N - 2)) 
  if (length(sel) > 0) inflMat[sel, 11] <- paste("*", inflMat[sel, 11], "*", sep = "")
 
  ## Crit 05: Cook's D
  COOK <- inflMat[, 12]
  sel <- which(COOK > qf(0.5, 2, N - 2))
  if (length(sel) > 0) inflMat[sel, 12] <- paste("*", inflMat[sel, 12], "*", sep = "")
  
  ## Crit 05: Leverage
  HAT <- inflMat[, 13]
  if (cutoff == "BKW") sel <- which(HAT > 2 * 2/N) else sel <- which(HAT > 3 * 2/N)
  if (length(sel) > 0) inflMat[sel, 13] <- paste("*", inflMat[sel, 13], "*", sep = "")
  
  ## Crit 06: Studentized residuals
  SR <- inflMat[, 15]
  sel <- which(abs(SR) > qt(0.975, N - 2 - 1))
  if (length(sel) > 0) inflMat[sel, 15] <- paste("*", inflMat[sel, 15], "*", sep = "")
 
  ## Crit 06: Hadi's measure
  HADI <- inflMat[, 14]
  sel <- which(HADI > median(HADI, na.rm = TRUE) + 2 * mad(HADI, na.rm = TRUE))
  if (length(sel) > 0) inflMat[sel, 14] <- paste("*", inflMat[sel, 14], "*", sep = "")
  
  ## Crit 07: CDR
  CDR <- inflMat[, 21]
  sel <- which(CDR > qbeta(0.05, 1, (N - 2 - 2)/2)/qbeta(0.05, 1, (N - 2 - 1)/2))
  if (length(sel) > 0) inflMat[sel, 21] <- paste("*", inflMat[sel, 21], "*", sep = "")
  
  ## Crit 08: Pena's Si
  SI <- inflMat[, 22]
  sel <- which(SI >= 0.9)
  if (length(sel) > 0) inflMat[sel, 22] <- paste("*", inflMat[sel, 22], "*", sep = "")
  
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