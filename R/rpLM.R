rpLM <- function(model, alpha = 0.05, R = 10000, plot = TRUE, verbose = TRUE, ...) {
  
  if (!is(model, "lm")) stop("'x' must be an univariate 'lm' object!") 
  pOrig <- summary(model)$coefficients[2, 4]
  
  if (verbose) cat("Running non-parametric cases bootstrap...\n")
  BOOT1 <- bootLM(model, type = "cases", alpha = alpha, R = R, ret.models = TRUE, ...) 
  if (verbose) cat("Running non-parametric residuals bootstrap...\n")
  BOOT2 <- bootLM(model, type = "residuals", alpha = alpha, R = R, ret.models = TRUE, ...) 
  if (verbose) cat("Running parametric residuals bootstrap...\n")
  BOOT3 <- bootLM(model, type = "parametric", alpha = alpha, R = R, ret.models = TRUE, ...) 
  
  pVec1 <- sapply(BOOT1, function(x) summary(x)$coef[2, 4])
  pVec2 <- sapply(BOOT2, function(x) summary(x)$coef[2, 4])
  pVec3 <- sapply(BOOT3, function(x) summary(x)$coef[2, 4])
  
  pAll <- data.frame("np.cases" = pVec1, "np.resid" = pVec2, "p.resid" = pVec3)
  
  if (plot) {
    stripchart(pAll, method = "jitter", vertical = TRUE, pch = 16, col = "#00000088", ylab = "P-value", log = "y", jitter = 0.3)
    abline(h = pOrig, lty = 1, lwd = 2, col = "darkgreen")
    abline(h = alpha, lty = 1, lwd = 2, col = "darkred")
    legend("bottomleft", c("initial p-value", "alpha"), fill = c("darkgreen", "darkred"))
  }
  
  if (pOrig <= alpha) SUM <- apply(pAll, 2, function(x) sum(x <= alpha, na.rm = TRUE)/R)
  else  SUM <- apply(pAll, 2, function(x) sum(x > alpha, na.rm = TRUE)/R)
  
  if (pOrig <= alpha & verbose) cat("\nOriginal P-value: ", signif(pOrig, 3), ".\nResults refer to proportion of bootstraps n(P < ", alpha, ").\n\n", sep = "")
  if (pOrig > alpha & verbose) cat("\nOriginal P-value: ", signif(pOrig, 3), ".\nResults refer to proportion of bootstraps n(P > ", alpha, ").\n\n", sep = "")
  return(round(SUM, 3))
}

