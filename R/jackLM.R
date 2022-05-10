## Jackknife function according to Quenouille 1956
jackLM <- function(formula, data = NULL, alpha = 0.05) {
  lmOrig <- lm(formula, data = data)
  DATA <- model.frame(lmOrig)
  N <- nrow(DATA)
  coefOrig <- coef(lmOrig)
  coefMat <- pseudoMat <- matrix(NA_real_, nrow = N, ncol = 2)
  
  for (i in 1:N) {
    lmJack <- lm(formula, data = DATA[-i, ])
    coefJack <- coef(lmJack)
    coefMat[i, ] <- coefJack
    pseudoMat[i, ] <- (N * coefOrig) - ((N - 1) * coefJack)
  }
  
  COEF <- apply(pseudoMat, 2, function(x) mean(x, na.rm = TRUE))
  SD <- apply(pseudoMat, 2, function(x) sd(x, na.rm = TRUE))
  SE <- SD/sqrt(N)
  PVAL <- 2 * pt(abs(COEF/SE), N - 1, lower.tail = FALSE)
  Lower <- COEF + SE * qt(alpha/2, df.residual(lmOrig))
  Upper <- COEF + SE * qt(1 - (alpha/2), df.residual(lmOrig))
  results <- data.frame('Estimate' = COEF, 'Std Error' = SE, `Lower bound` = Lower, `Upper bound` = Upper, `p-value` = PVAL)
  rownames(results) <- names(coefOrig)
  return(results)
}
