lmExact <- function(
  x = 1:20, 
  ny = 1,
  intercept = 0, 
  slope = 0.1, 
  error = 0.1, 
  seed = 123,
  pval = NULL,
  rsq = NULL,
  plot = TRUE, 
  verbose = FALSE,
  ...) 
{
  ## create linear model environment
  lmEnv <- new.env()
  
  ## define x-vector
  x <- rep(x, each = ny)
  
  ## create error vector, either s.d. fraction or own vector
  if (length(error) == 1) {
    set.seed(seed)
    errorVec <- rnorm(length(x), 0, error)
  } else {
    error <- as.numeric(error)
    if (length(error) != length(x)) stop("'error' must have same length as 'x'!")
    errorVec <- error
  }
  
  ## define optimizing function for slope
  optFct <- function(slope) {
    y <- intercept + slope * x + errorVec
    LM <- lm(y ~ x, ...)
    RESID <- residuals(LM)
    y <- intercept + slope * x + RESID
    LM <- lm(y ~ x, ...)
    PVAL <- summary(LM)$coefficients[2, 4]
    OUT <- abs(PVAL - pval)
    assign("LM", LM, envir = lmEnv)
    return(OUT)
  }
  
  ## if p-value is given, optimize slope, else return intercept + slope * x + residuals 
  if (!is.null(pval)) {
    OPT <- suppressWarnings(optim(slope, optFct, control = list(trace = 0, maxit = 500)))
    LM <- get("LM", envir = lmEnv)
    if(is.character(all.equal(pval, summary(LM)$coefficients[2, 4], tolerance = 0.01 * pval))) 
      print("Warning: Optimizer has not converged to the desired p-value! Try a different 'slope' starting value.")
  } else
    if (!is.null(rsq)) {
      ## This one is taken from http://stats.stackexchange.com/questions/15011/
      ## generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
      if (rsq < 0 | rsq > 1) stop("R-square must be in [0, 1] !")
      y <- intercept + slope * x + errorVec
      rho   <- sqrt(rsq)        
      theta <- acos(rho)  
      MAT     <- cbind(x, y)  
      scaleMAT  <- scale(MAT, center = TRUE, scale = FALSE)
      Id   <- diag(length(x))                               
      Q    <- qr.Q(qr(scaleMAT[ , 1, drop=FALSE]))  
      P    <- tcrossprod(Q)      
      x2o  <- (Id-P) %*% scaleMAT[ , 2]    
      Xc2  <- cbind(scaleMAT[ , 1], x2o)   
      tempY    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2))) 
      y <- tempY[ , 2] + (1 / tan(theta)) * tempY[ , 1] 
      LM <- lm(y ~ x, ...)
      INTER <- coef(LM)[1]
      y <- y + (intercept - coef(LM)[1])
      LM <- lm(y ~ x, ...)
      assign("LM", LM, envir = lmEnv)
    }
  else {
    y <- intercept + slope * x + errorVec
    LM <- lm(y ~ x, ...)
    RESID <- residuals(LM)
    y <- intercept + slope * x + RESID
    LM <- lm(y ~ x, ...)
    assign("LM", LM, envir = lmEnv)
  }
  
  LM <- get("LM", envir = lmEnv)
  
  ## optional plotting
  if (plot) {
    DATA <- LM$model
    plot(DATA[, 2], DATA[, 1], pch = 16, xlab = "x", ylab = "y", las = 1, ...)
    abline(LM, ...)
    grid()
  }
  
  ## display summary
  if (verbose) print(summary(LM))
  
  return(list(lm = LM, x = LM$model[, 2], y = LM$model[, 1], summary = summary(LM)))
}