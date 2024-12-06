pcomp <- function(x, y = NULL, R = 1000, alpha = 0.05, ...) {
  
  ## use either model or raw data
  if (is(x, "lm")) {
    xorig <- x
    x <- model.frame(xorig)[, 2]
    y <- model.frame(xorig)[, 1]
  }
  
  if (!is(x, "lm") & is.null(y)) stop("If 'x' is not an 'lm' model, 'x' and 'y' raw data must be supplied!")
  
  DATA <- data.frame(x = x, y = y)
  
  ## Original model
  LM1 <- lm(y ~ x, data = DATA, ...)
  p1 <- summary(LM1)$coefficients[2, 4]
  
  ## Leave strongest outlier out
  LM2 <- lmInfl(LM1, verbose = FALSE, ...)
  if (LM2$origP <= alpha) p2 <- max(LM2$raw[, "looP"], na.rm = TRUE) else p2 <- min(LM2$raw[, "looP"], na.rm = TRUE)
  
  ## MM estimator
  LM3 <- lmrob(y ~ x, data = DATA, control = lmrob.control(max.it = 500), ...)
  p3 <- summary(LM3)$coefficients[2, 4]
  
  ## Theil-Sen estimator
  p4 <- cor.test(x, y, method = "kendall", exact = TRUE, continuity = TRUE, ...)$p.value
  
  ## Least absolute deviations regression
  LM5 <- lad(y ~ x, ...)
  p5 <- summary(LM5)$coefficients[2, 4]
  
  ## quantile regression
  LM6 <- rq(y ~ x, tau = 0.9, ...)
  p6 <- summary(LM6, se = "ker")$coefficients[2, 4]
  
  ## weighted regression with isotree 1/scores as weights
  IF <- isolation.forest(matrix(residuals(LM1)), ndim = 1, ntrees = 1000, nthreads = 1, ...)
  scores <- predict(IF, matrix(residuals(LM1)), type = "score")
  LM7 <- lm(y ~ x, weights = 1/scores)
  p7 <- summary(LM7)$coefficients[2, 4]
  
  ## Bootstrap
  LM8 <- bootLM(LM1, R = R, ...)
  p8 <- LM8[2, 5]
  
  ## Jackknife
  LM9 <- jackLM(LM1, ...)
  p9 <- LM9[2, 5]
  
  return(c("p.lm" = p1, "p.infl" = p2, "p.rob" = p3, "p.ts" = p4, "p.lad" = p5, "p.quant" = p6, "p.if" = p7, "p.boot" = p8, "p.jack" = p9))
}
