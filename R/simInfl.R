simInfl <- function(
  x = 1:10, 
  slope = 0.02, 
  intercept = 1,
  error = 0.05, 
  nrev = 1000, 
  ...)
{
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
  
  ## preallocate result matrix
  MAT <- matrix(NA, nrow = 1000000, ncol = 20)
  
  ## initialize reversal counter
  isRev <- 0
  
  ## preallocate result vectors
  seedVec <- pVec <- rep(NA, 1000000)
  mList <- vector("list", length = 1000000)
  
  ## loop until isRev = nrev
  for (i in 1:1000000) {
    ## create exact model 
    LME <- lmExact(x = x, slope = slope, intercept = intercept, 
                   error = error, plot = FALSE, seed = i, verbose = FALSE, ...)
    
    ## get reversal results
    RES <- lmInfl(LME$lm, verbose = FALSE, ...)
    
    ## populate result matrix and increase counters
    seedVec[i] <- i
    pVec[i] <- RES$origP
    if (!is.null(RES$sel)) {
      MAT[i, ] <- as.numeric(RES$infl[RES$sel[1], ])
      isRev <- isRev + 1
      counter(isRev)
    }
    
    ## store all models
    mList[[i]] <- LME$lm
    
    ## break loop if 'nrev' reversals are counted
    if (isRev == nrev) {
      niter <- i
      break
    }
  }
 
  ## create result list
  colnames(MAT) <- colnames(RES$infl)
  MAT <- cbind(seed = seedVec, origP = pVec, MAT)
  
  ## condense result matrix & model list
  sel <- which(is.na(MAT[, "Idx"]))
  MAT <- MAT[-sel, ]
  mList <- mList[1:i]
  
  OUT <- list(models = mList, mat = MAT)
  return(OUT)
}