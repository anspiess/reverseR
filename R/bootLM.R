bootLM <- function (formula, data = NULL, R = 10000, alpha = 0.05)  
{
  DATA <- data
  LM <- lm(formula, data = DATA)
  DATA <- model.frame(LM)
  
  bf <- function(formula, data, indices) {
    d <- data[indices, ]
    model <- lm(formula, data = d)
    return(coef(model))
  }
  
  b <- boot(data = DATA, statistic = bf, R = R, formula = formula, parallel = "no", ncpus = 1, cl = NULL)
  estimate <- apply(b$t, 2L, mean, na.rm = TRUE)
  stderr <- apply(b$t, 2L, sd, na.rm = TRUE)
  p <- length(coef(LM))
  results <- data.frame('Estimate' = estimate, 'Std Error' = stderr, `Lower bound` = rep(NA, p), `Upper bound` = rep(NA, p), `p-value` = rep(NA, p))
  for (i in 1:p) {
    ci <- boot.ci(b, conf = 1 - alpha, type = "perc", index = i)
    results[i, 3:4] <- ci$percent[, 4:5]
  }
  for (i in 1:p) {
    results[i, 5] <- boot.pval(b, type = "perc", theta_null = 0, pval_precision = 1/R, index = i)
  }
  rownames(results) <- names(b$t0)
  return(results)
}
