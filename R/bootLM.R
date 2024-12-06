bootLM <- function(model, type = c("cases", "residuals", "residuals2", "parametric"), R = 10000, alpha = 0.05, ret.models = FALSE)  
{
  type <- match.arg(type)
  lmOrig <- model 
  DATA <- model.frame(lmOrig)
  FITTED <- fitted(lmOrig)
  RESID <- residuals(lmOrig)
  FORMULA <- formula(lmOrig)
  SIGMA <- summary(lmOrig)$sigma
  set.seed(NULL)
  
  if (type == "cases") {
    modList <- vector("list", length = R)
    for (i in 1:R) {
      sel <- sample(1:nrow(DATA), replace = TRUE)
      d <- DATA[sel, ]
      modList[[i]] <- lm(FORMULA, data = d)
    }
  }
  
  if (type == "residuals") {
    modList <- vector("list", length = R)
    for (i in 1:R) {
      sel <- sample(1:nrow(DATA), replace = TRUE)
      d <- DATA; d[, 1] <- FITTED + RESID[sel]
      modList[[i]] <- lm(FORMULA, data = d)
    }
  }
  
  if (type == "residuals2") {
    modList <- vector("list", length = R)
    RESID2 <- RESID/sqrt(1 - hatvalues(lmOrig)) - mean(RESID, na.rm = TRUE)
    for (i in 1:R) {
      sel <- sample(1:nrow(DATA), replace = TRUE)
      d <- DATA; d[, 1] <- FITTED + RESID2[sel]
      modList[[i]] <- lm(FORMULA, data = d)
    }
  }
  
  if (type == "parametric") {
    modList <- vector("list", length = R)
    for (i in 1:R) {
      d <- DATA; d[, 1] <- FITTED + rnorm(nrow(DATA), 0, SIGMA)
      modList[[i]] <- lm(FORMULA, data = d)
    }
  }
  
  COEF <- t(sapply(modList, coef))
  estimate <- apply(COEF, 2L, mean, na.rm = TRUE)
  stderr <- apply(COEF, 2L, sd, na.rm = TRUE)
  p <- length(coef(lmOrig))
  results <- data.frame('Estimate' = estimate, 'Std Error' = stderr, `Lower bound` = rep(NA, p), `Upper bound` = rep(NA, p), `p-value` = rep(NA, p))
  b <- list(t0 = coef(lmOrig), t = COEF, R = R, data = DATA)
  
  for (i in 1:p) {
    ci <- boot.ci(b, conf = 1 - alpha, type = "perc", index = i)
    results[i, 3:4] <- ci$perc[, 4:5]
  }
  for (i in 1:p) {
    results[i, 5] <- boot.pval(b, type = "perc", theta_null = 0, pval_precision = 1/R, index = i)
  }
  
  rownames(results) <- names(b$t0)
  if (!ret.models) return(results) else return(modList)
}
