\name{lmThresh}
\alias{lmThresh}
\encoding{latin1}

\title{Finds and analyzes significance reversal regions for each response value}

\description{
This function finds (by iterating through a grid of values for each response) the approximate response value range(s) in which the regression is significant (when inside) or not (when outside), as defined by \code{alpha}. Here, two scenarios can be tested: i) if \code{newobs = FALSE} (default), the model's significance is tested by shifting \eqn{y_i} along the search grid. If \code{newobs = TRUE}, \eqn{y_i} is kept fixed and a \code{new} \code{obs}ervation \eqn{y_{2i}} is added and shifted along the search grid. Hence, this function tests the regression for the sensitivity of being reversed in its significance through minor shifting of the original or added response values, as opposed to the effect of point removal (\code{\link{lmInfl}}).
}

\usage{
lmThresh(model, factor = 5, alpha = 0.05, 
         method = c("pearson", "spearman"),
         steps = 10000, newobs = FALSE, ...) 
}

\arguments{
\item{model}{the linear model of class \code{\link{lm}}.}
\item{factor}{a factor for the initial search grid. See 'Details'.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{method}{select either parametric (\code{"pearson"}) or rank-based (\code{"spearman"}) statistics.}
\item{steps}{the number of steps within the search range. See 'Details'.}
\item{newobs}{logical. Should the significance region for each \eqn{y_i} be calculated from shifting \eqn{y_i} or from keeping \eqn{y_i} fixed and adding a new observation \eqn{y2_i}?}
\item{...}{other arguments to future methods.}
}

\details{
In a first step, a grid is created with a range from \eqn{y_i \pm \mathrm{factor} \cdot \mathrm{range}(y_{1...n})} with \code{steps} cuts. For each cut, the \emph{p}-value is calculated for the model when \eqn{y_i} is shifted to that value (\code{newobs = FALSE}) or a second observation \eqn{y_{2i}} is added to the fixed \eqn{y_i} (\code{newobs = TRUE}). When the original model \eqn{y =  \beta_0 + \beta_1x + \varepsilon} is significant (\emph{p} < \code{alpha}), there are two boundaries that result in insignificance: one decreases the slope \eqn{\beta_1} and the other inflates the standard error \eqn{\mathrm{s.e.}(\beta_1)} in a way that \eqn{P_t(\frac{\beta_1}{\mathrm{s.e.}(\beta_1)}, n-2) > \alpha}. If the original model was insignificant, also two boundaries exists that either increase \eqn{\beta_1} or reduce \eqn{\mathrm{s.e.}(\beta_1)}. Often, no boundaries are found and increasing the \code{factor} grid range may alleviate this problem.

This function is quite fast (~ 300ms/10 response values), as the slope's \emph{p}-value is calculated from the \code{corr.test} function of the 'psych' package, which utilizes matrix multiplication and vectorized \code{\link{pt}} calculation. The vector of correlation coefficients \eqn{r_i} from the \code{\link{cor}} function is transformed to t-values by \deqn{t_i = \frac{r_i\sqrt{n-2}}{\sqrt{1-r_i^2}}} which is equivalent to that employed in the linear regression's slope test. 
}

\value{
A list with the following items:\cr
\item{x}{the predictor values.}
\item{y}{the response values.}
\item{pmat}{the \emph{p}-value matrix, with \code{length(x)} columns and \code{steps} rows.}
\item{alpha}{the selected \eqn{\alpha}-level.}
\item{ySeq}{the grid sequence for which the algorithm calculates \emph{p}-values when \eqn{y_i} is shifted within.}
\item{model}{the original \code{\link{lm}} model.}
\item{data}{the original \code{\link{model.frame}}.}
\item{eosr}{the y-values of the ends of the significance region.}
\item{diff}{the \eqn{\Delta} value between \eqn{y_i} and the nearest border of significance reversal.}
\item{closest}{the (approx.) value of the nearest border of significance reversal.}
\item{newobs}{should a new observation be added?}
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Significant model, no new observation.
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(length(a), 0, 1)
LM1 <- lm(b ~ a)
res1 <- lmThresh(LM1)
threshPlot(res1)
stability(res1)

## Insignificant model, no new observation.
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(length(a), 0, 2)
LM2 <- lm(b ~ a)
res2 <- lmThresh(LM2)
threshPlot(res2)
stability(res2)

## Significant model, new observation.
## Some significance reversal regions
## are within the prediction interval,
## e.g. 1 to 6 and 14 to 20.
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(length(a), 0, 1)
LM3 <- lm(b ~ a)
res3 <- lmThresh(LM3, newobs = TRUE)
threshPlot(res3)
stability(res3)

## More detailed example to the above:
## a (putative) new observation within the
## prediction interval may reverse significance.
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(length(a), 0, 1)
LM1 <- lm(b ~ a)
summary(LM1) # => p-value = 0.02688
res1 <- lmThresh(LM1, newobs = TRUE)
threshPlot(res1)
st <- stability(res1, pval = TRUE)
st$stats # => upper prediction boundary = 7.48
         # and eosr = 6.49
stabPlot(st, 1)
## reverse significance if we add a new response y_1 = 7
a <- c(1, a)
b <- c(7, b)
LM2 <- lm(b ~ a)
summary(LM2) # => p-value = 0.0767
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

