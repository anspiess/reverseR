lmThresh <- function(model, factor = 5, alpha = 0.05, method = c("pearson", "spearman"), steps = 10000, newobs = FALSE, ...) {
  ## define correlation function for fast p-value calculation
  ## taken from the 'corr.test' function of the 'psych' package
  corr.test <- function (x, y) 
  {
    if (method == "pearson") r <- cor(x, y, use = "complete.obs", method = "pearson")
    else  r <- cor(x, y, use = "complete.obs", method = "spearman")
    n <- t(!is.na(x)) %*% (!is.na(y)); n <- min(n)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- 2 * (1 - pt(abs(t), (n - 2)))
    se <- sqrt((1 - r * r)/(n - 2))
    return(list(p = p, t = t, se = se, r = r))
  }
  
  method <- match.arg(method)
  
  ## counter function
  counter <- function (i) 
  {
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
  
  if (class(model) != "lm") stop("Not an 'lm' model!")
  
  ## p-value of full model
  fullModelPval <- summary(model)$coefficients[2, 4]
  
  ## extract data and define range
  DATA <- eval(model$model)
  X <- DATA[, 2]; Y <- DATA[, 1]
  MAX <- max(Y, na.rm = TRUE); MIN <- min(Y, na.rm = TRUE)
  RANGE <- MAX-MIN
  
  ## create grid sequence
  SEQ <- seq(MIN - factor * RANGE, MAX + factor * RANGE, length.out = steps)
  
  ## preallocate result vector and list
  pMat <- matrix(NA, ncol = nrow(DATA), nrow = steps)
  
  ## create replicated y-matrix
  yMat <- replicate(steps, Y)
  
  cat("Calculating p-values within [", min(SEQ, na.rm = TRUE), ";", max(SEQ, na.rm = TRUE), "]\n", sep = "")
  ## iterate over x-range and find y-value that delivers 'alpha' p-value
  for (i in 1:nrow(DATA)) {
    
    ## reinitialize y-matrix
    tempMat <- yMat
    
    ## initialize X2
    X2 <- X
    
    ## fill with sequence at y_i
    tempMat[i, ] <- SEQ
    
    ## if 'newobs' = TRUE, add original y_i value
    if (newobs) { 
      tempMat <- rbind(rep(Y[i], steps), tempMat)
      X2 <- c(X2[i], X2)
    } 
   
    ## calculate p-values for grid sequence
    counter(i)
    pVec <- as.numeric(corr.test(X2, tempMat)$p)
    
    ## put a p = 1 at both ends to create bordered region
    pVec[1] <- pVec[length(pVec)] <- 1
    
    ## insert into matrix
    pMat[, i] <- rev(pVec)
  }
  cat("\n")
 
  ## row/col names for pMat
  RN <- rownames(pMat) <- rSEQ <- rev(SEQ); colnames(pMat) <- 1:length(X)
  
  ## preallocate vectors and lists
  vecDIFF <- vecCLOSE <- rep(NA, length(Y))
  IntMat <- matrix(NA, ncol = 2, nrow = length(Y))
  
  ## calculate delta of y_i and end of significance region
  for (i in 1:ncol(pMat)) {
    logVec <- rev(pMat[, i] <= alpha) ## get significance region
    diffVec <- diff(logVec) ## calculate delta's
    selVec <- which(diffVec == 1 | diffVec == -1) ## get border indices of significance region
    selVec[1] <- selVec[1] + 1 ## need to increment by 1
    if (length(selVec) > 2) selVec[3] <- selVec[3] + 1  
    Int <- rev(RN[selVec]) ## get rownames = y-values of borders
    Int <- as.numeric(names(selVec))
    
    ## decide, if two significance are present
    if (length(Int) == 4) {
      if (Y[i] > Int[1] & Y[i] < Int[2]) Int <- Int[1:2]
      if (Y[i] > Int[3] & Y[i] < Int[4]) Int <- Int[3:4]
    }
    
    IntMat[i, ] <- Int
    selInt <- findInterval(Y[i], Int) ## in which region is y_i located?
    
    ## Case 1: y_i < lower
    if (selInt == 0) {
      vecDIFF[i] <- Int[1] - Y[i]
      vecCLOSE[i] <- Int[1]
    }
    
    ## Case 2: y_i > upper
    else if (length(Int) == 2 & selInt == 2) {
      vecDIFF[i] <- Y[i] - Int[2] 
      vecCLOSE[i] <- Int[2]
    }
    
    ##  Case 3: lower < y_i < upper
    else {
      vecDIFF[i] <- min(abs(Y[i] - c(Int[selInt], Int[selInt + 1])))
      vecCLOSE[i] <- Int[which.min(abs(Y[i] - c(Int[selInt], Int[selInt + 1])))]
    }
  }
  
  OUT <- list(x = X, y = Y, pmat = pMat, alpha = alpha, ySeq = SEQ, model = model,  
              data = DATA, eosr = IntMat, diff = vecDIFF, closest = vecCLOSE, newobs = newobs)
  class(OUT) <- "threshsearch"
  return(OUT)
}

threshPlot <- function(thresh, bands = FALSE, ...) {
  ## check for correct class
  if (class(thresh) != "threshsearch") stop("object is not a result of 'lmThresh' !")
  
  ## extract x/y, p-value matrix, sequence, alpha, model
  X <- thresh$x; Y <- thresh$y; pMat <- thresh$pmat; SEQ <- thresh$ySeq; 
  alpha <- thresh$alpha; Model <- thresh$model
  
  ### calculate full model p-value
  fullModelPval <- summary(thresh$model)$coefficients[2, 4]
  
  ## plot x/y in 'lmThresh's range and regression line 
  plot(X, Y, xlab = names(thresh$st$stats)[1], ylab = names(thresh$st$stats)[2], 
       ylim = c(min(SEQ), max(SEQ)), type = "n", ...)
  grid()
  
  ## add significance regions
  arrows(X, thresh$eosr[, 1], X, thresh$eosr[, 2], length = 0, col = "green4", lwd = 3)
  
  ## add regression line and points
  abline(thresh$model, lwd = 2, col = "darkred", ...)
  points(X, Y, pch = 16, cex = 1.2, ...)
  
  ## add confidence and prediction line
  if (bands | thresh$newobs) {
    aT <- par()$usr
    xSeq <- seq(aT[1], aT[2], length.out = 100)
    xDat <- data.frame(xSeq); colnames(xDat) <- colnames(Model$model)[2]
    CONF <- predict(Model, interval = "confidence", level = 1 - alpha, ...)
    lines(X, CONF[, 2], lty = 1, lwd = 1.5, type = "s", ...)
    lines(X, CONF[, 3], lty = 1, lwd = 1.5, type = "s", ...)
    PRED <- predict(Model, interval = "prediction", level = 1 - alpha, ...)
    lines(X, PRED[, 2], lty = 2, lwd = 1.5, type = "s", ...)
    lines(X, PRED[, 3], lty = 2, lwd = 1.5, type = "s", ...)
  }
}
