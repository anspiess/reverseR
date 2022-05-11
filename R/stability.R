stability <- function(x, pval = FALSE, ...) {
  
  ## for results of 'lmInfl'
  if (as.character(class(x)) == "influencers") {
    stab <- 1 - (length(x$sel)/nrow(model.frame(x$origModel)))
    logVec <- rep(F, nrow(model.frame(x$origModel)))
    logVec[x$sel] <- TRUE
    statDat <- NULL
  }
  
  ## for results of 'lmThresh'
  else if (as.character(class(x)) == "threshsearch") {
    ## get original data
    Model <- x$model
    DATA <- model.frame(Model); X <- DATA[, 2]; Y <- DATA[, 1]
    eosr <- x$eosr
    alpha <- x$alpha
    N <- length(X)
    DF <- df.residual(Model)
    Fitted <- fitted(Model)
    
    ## get prediction interval
    Pred <- suppressWarnings(predict(Model, interval = "prediction", level = 1 - alpha, ...))
   
    ## check if closest (closest end of significance region) is within prediction interval
    isWithin1 <- ifelse(x$eosr[, 1] > Pred[, 2] & x$eosr[, 1] < Pred[, 3], TRUE, FALSE)
    isWithin2 <- ifelse(x$eosr[, 2] > Pred[, 2] & x$eosr[, 2] < Pred[, 3], TRUE, FALSE)
    stab <- 1 - sum(isWithin1 | isWithin2, na.rm = TRUE)/length(Y)
    logVec <- isWithin1 | isWithin2
    
    ## if p-value, calculate the integral from the closest point 
    ## to the nearest prediction interval border
    if (pval) {
      ## define scaled / shifted density function of t-distribution
      dt.scaled <- function(x, df, mu, s) 1/s * dt((x - mu)/s, df)
      
      ## calculate "standard error of prediction
      MSE <- sum(residuals(Model)^2)/DF
      SE <- sqrt(MSE * (1 + 1/N + (X - mean(X))^2/sum((X - mean(X))^2)))
      
      ## calculate upper/lower prediction interval
      Upper <- Fitted + qt(1-alpha/2, DF) * SE
      Lower <- Fitted - qt(1-alpha/2, DF) * SE
      
      ## preallocate probability vector
      pval <- rep(0, length(Y))
      
      ## get probabilities as the integral between "eosr" (significance region end) and
      ## the 1-alpha/2 prediction interval upper/lower bound
      for (i in 1:length(Y)) {
        if (isWithin1[i]) pval[i] <- integrate(dt.scaled, lower = Lower[i], upper = x$eosr[i, 1], 
                                            df = DF, mu = Fitted[i], s = SE[i])$value
        if (isWithin2[i]) pval[i] <- integrate(dt.scaled, lower = x$eosr[i, 2], upper = Upper[i], 
                                                     df = DF, mu = Fitted[i], s = SE[i])$value
      }
      statDat <- data.frame(x = X, y = Y, fitted = Fitted, eosr = x$eosr, 
                            lower = Lower, upper = Upper, se = SE, pval = pval)
    } else statDat <- NULL
  }
  
  cat("Stability Index:", stab, "\n\n")
  invisible(list(stab = stab, within = logVec, stats = statDat))
}

stabPlot <- function(stab, which = NULL, ...) {
  ## define scaled/shifted t-distribution
  dt.scaled <- function(x, df, mu, s) 1/s * dt((x - mu)/s, df)
  
  ## check type
  if (is.null(stab$stats)) stop("This is not a result from 'stability(lmThresh-object, pval = TRUE)'!")
  # check response value selection
  if (is.null(which)) stop("Please select a response value to plot.")
  ## check proper selection
  stats <- stab$stats
  which <- which[1]
  if (!(which %in% 1:nrow(stats))) stop("'which' must be within 1...", nrow(stats), "!")
  
  Sel <- stats[which, ]
  
  ## plot density curve
  usr <- par()$usr
  MIN <- Sel[, "fitted"] - 1.5 * (Sel[, "fitted"] - Sel[, "lower"])
  MAX <- Sel[, "fitted"] + 1.5 * (Sel[, "upper"] - Sel[, "fitted"])
  SEQ <- seq(from = MIN, to = MAX, length.out = 1000)
  DENS <- dt.scaled(SEQ, nrow(stats) - 2, Sel[, "fitted"], Sel[, "se"])
  plot(SEQ, DENS, type = "l", lwd = 2, xlab = "Y value", ylab = "Density")
  
  ## fill 1-alpha tails
  SEQ1 <- seq(from = MIN, to = Sel[, "lower"], length.out = 1000)
  DENS1 <- dt.scaled(SEQ1, nrow(stats) - 2, Sel[, "fitted"], Sel[, "se"])
  segments(SEQ1, 0, SEQ1, DENS1, col = "darkred", lwd = 2)
  SEQ2 <- seq(from = Sel[, "upper"], to = MAX, length.out = 1000)
  DENS2 <- dt.scaled(SEQ2, nrow(stats) - 2, Sel[, "fitted"], Sel[, "se"])
  segments(SEQ2, 0, SEQ2, DENS2, col = "darkred", lwd = 2)
  
  ## include significance region from lmThresh
  rect(Sel[, "eosr.1"], 0, Sel[, "eosr.2"], 0.003, col = "green4")
  
  ## include y-value
  points(Sel[, "y"], 0, pch = 16, cex = 2, xpd = TRUE)
  
  ## plot delta-region
  if (Sel[, "eosr.2"] > Sel[, "lower"] & Sel[, "eosr.2"] < Sel[, "upper"]) {
    SEQ3 <- seq(from = Sel[, "eosr.2"], to = Sel[, "upper"], length.out = 1000)
    DENS3 <- dt.scaled(SEQ3, nrow(stats) - 2, Sel[, "fitted"], Sel[, "se"])
    segments(SEQ3, 0, SEQ3, DENS3, col = "dodgerblue4", lwd = 2)
  }
  if (Sel[, "eosr.1"] > Sel[, "lower"] & Sel[, "eosr.1"] < Sel[, "upper"]) {
    SEQ3 <- seq(from = Sel[, "lower"], to = Sel[, "eosr.1"], length.out = 1000)
    DENS3 <- dt.scaled(SEQ3, nrow(stats) - 2, Sel[, "fitted"], Sel[, "se"])
    segments(SEQ3, 0, SEQ3, DENS3, col = "dodgerblue4", lwd = 2)
  }
  
  ## add stats in title
  STR <- paste("Y:", round(Sel[, "y"], 2), "    Fitted:", round(Sel[, "fitted"], 2),
               "    EOSR: [", round(Sel[, "eosr.1"], 2), "/", round(Sel[, "eosr.2"], 2), "]\n",
               "Lower:", round(Sel[, "lower"], 2), "    Upper", round(Sel[, "upper"], 2),
               "    P-value:", signif(Sel[, "pval"], 3))
  title(main = STR)
}

