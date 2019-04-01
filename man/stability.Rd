\name{stability}
\alias{stability}
\encoding{latin1}

\title{Calculates stability values for results of 'lmInfl', 'lmMult' and 'lmThresh'}

\description{
This function calculates stability values for LOO (\code{\link{lmInfl}}), LMO (\code{\link{lmMult}}) and response value shifting/addition (\code{\link{lmThresh}}).
}

\usage{
stability(x, pval = FALSE, ...) 
}

\arguments{
\item{x}{a result of either \code{\link{lmInfl}}, \code{\link{lmMult}} or \code{\link{lmThresh}}.}
\item{pval}{logical. If \code{TRUE}, for \code{\link{lmThresh}}, objects an exact \emph{p}-value is calculated for a future response to reverse significance.}
\item{...}{other parameters, not yet implemented.}
}

\details{
For results of \code{\link{lmInfl}}:\cr
A [0, 1]-bounded stability measure \eqn{S = 1-\frac{n}{N}}, with \eqn{n} = number of influencers (significance reversers) and \eqn{N} = total number of response values.\cr\cr
For results of \code{\link{lmMult}}:\cr
For each 1...\code{max}, the percentage of all resamples that did *NOT* result in significance reversal.\cr\cr
For results of \code{\link{lmThresh}}:\cr
A [0, 1]-bounded stability measure \eqn{S = 1-\frac{n}{N}}, with \eqn{n} = number of response values where one of the ends of the significance region is within the prediction interval and \eqn{N} = total number of response values.\cr
If \code{pval = TRUE}, the exact \emph{p}-value is calculated in the following manner:\cr

1) Mean square error (MSE) and prediction standard error (se) are calculated from the linear model:\cr
\deqn{\mathrm{MSE} = \sum_{i=1}^n \frac{(y_i - \hat{y}_i)^2}{n-2} \quad\quad \mathrm{se}_i = \sqrt{\mathrm{MSE} \cdot \left(1 + \frac{1}{n} + \frac{(x_i - \bar{x}_i)^2}{\sum_{i=1}^n (x_i - \bar{x}_i)^2}\right)}}
2) Upper and lower prediction intervals boundaries are calculated for each \eqn{\hat{y}_i}:\cr
\deqn{\hat{y}_i \pm Q_t(\alpha/2, n-2) \cdot \rm{se}_i}
The prediction interval around \eqn{\hat{y}_i} is a scaled/shifted \eqn{t}-distribution with density function \deqn{P_{tss}(y, n-2) = \frac{1}{\rm{se}_i} \cdot P_t\left(\frac{y - \hat{y}_i}{\rm{se}_i}, n-2\right)}, where \eqn{P_t} is the density function of the central, unit-variance \eqn{t}-distribution.\cr
3) The probability of either shifting the response value (if \code{lmThresh(..., newobs = FALSE)}) or including a future response value \eqn{y_{2i}} (if \code{lmThresh(..., newobs = TRUE)}) to reverse the significance of the linear model is calculated as the integral between the end of the significance region (eosr) and the upper/lower \eqn{\alpha/2, 1-\alpha/2} prediction interval:
\deqn{P(\mathrm{reverse}) = \int_{\mathrm{eosr}}^{1-\alpha/2} P_{tss}(y, n-2)dy \quad \mathrm{or} \quad \int_{\alpha/2}^{\mathrm{eosr}} P_{tss}(y, n-2)dy}
}

\value{
The stability value.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## See examples in 'lmInfl' and 'lmThresh'.

## The implemented strategy of calculating the
## probability of significance reversal, as explained above
## and compared to 'stabPlot'.
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(length(a), 0, 1)
LM1 <- lm(b ~ a)
res1 <- lmThresh(LM1, newobs = TRUE)
st1 <- stability(res1, pval = TRUE)

## Let's check that the prediction interval encompasses 95\%:
dt.scaled <- function(x, df, mu, s) 1/s * dt((x - mu)/s, df)
integrate(dt.scaled, lower = st1$stats[1, "lower"], st1$stats[1, "upper"], 
          df = 18, mu = st1$stats[1, "fitted"], s = st1$stats[1, "se"])
## => 0.95 with absolute error < 8.4e-09

## This is the interval between "end of significance region" and upper 
## prediction boundary:
integrate(dt.scaled, lower = st1$stats[1, "eosr.2"], st1$stats[1, "upper"], 
          df = 18, mu = st1$stats[1, "fitted"], s = st1$stats[1, "se"])
## => 0.09264124 with absolute error < 1e-15

## We can recheck this value by P(B) - P(A):
pt.scaled <- function(x, df, mu, s) pt((x - mu)/s, df)
pA <- pt.scaled(x = st1$stats[1, "eosr.2"], df =  18, mu = st1$stats[1, "fitted"], 
                s = st1$stats[1, "se"])
0.975 - pA 
##  => 0.09264124 as above
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

