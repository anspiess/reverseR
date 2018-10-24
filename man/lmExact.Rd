\name{lmExact}
\alias{lmExact}
\encoding{latin1}

\title{Create random values that deliver linear regressions with exact parameters}

\description{
Takes self-supplied or generated random values and transforms these as to deliver linear regressions \eqn{y = \beta_0 + \beta_1x + \varepsilon} (with potential replicates) with either\cr

\bold{1)} exact slope \eqn{\beta_1} and intercept \eqn{\beta_0},\cr
\bold{2)} exact \emph{p}-value and intercept \eqn{\beta_0}, or\cr
\bold{3)} exact \eqn{R^2} and intercept \eqn{\beta_0}.\cr

Intended for testing and education, not for cheating ! ;-)
}

\usage{
lmExact(x = 1:20, ny = 1, intercept = 0, slope = 0.1, error = 0.1, 
        seed = 123, pval = NULL, rsq = NULL, plot = TRUE, verbose = FALSE, ...) 
}

\arguments{
\item{x}{the predictor values.}
\item{ny}{the number of replicate response values per predictor value.}
\item{intercept}{the desired intercept \eqn{\beta_0}.}
\item{slope}{the desired slope \eqn{\beta_1}.}
\item{error}{if a single value, the standard deviation \eqn{\sigma} for sampling from a normal distribution, or a user-supplied vector of length \code{x} with random deviates.}
\item{seed}{the random generator seed for reproducibility.}
\item{pval}{the desired \emph{p}-value of the slope.}
\item{rsq}{the desired \eqn{R^2}.}
\item{plot}{logical. If \code{TRUE}, the linear regression is plotted.}
\item{verbose}{logical. If \code{TRUE}, a summary is printed to the console.}
\item{...}{other arguments to \code{\link{lm}} or \code{\link{plot}}.}
}

\details{
For case \bold{1)}, the \code{error} values are added to the exact \eqn{(x_i, \beta_0 + \beta_1 x_i)} values, the linear model \eqn{y_i = \beta_0 + \beta_1 x_i + \varepsilon} is fit, and the residuals \eqn{y_i - \hat{y_i}} are re-added to \eqn{(x_i, \beta_0 + \beta_1 x_i)}.\cr
For case \bold{2)}, the same as in \bold{1)} is conducted, however the slope delivering the desired \emph{p}-value is found by an optimizing algorithm.\cr
Finally, for case \bold{3)}, a QR reconstruction, rescaling and refitting is conducted, using the code found under 'References'.
}

\value{
A list with the following items:\cr
\item{lm}{the linear model of class \code{lm}.} 
\item{x}{the predictor values.}
\item{y}{the (random) response values.}
\item{summary}{the model summary for quick checking of obtained parameters.}

Using both \code{x} and \code{y} will give a linear regression with the desired parameter values when refitted.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
For method \bold{3}):\cr
http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable.
}

\examples{
## No replicates, intercept = 3, slope = 0.2, sigma = 2, n = 20.
res1 <- lmExact(x = 1:20, ny = 1, intercept = 3, slope = 2, error = 2)

## Same as above, but with 3 replicates, sigma = 1,  n = 20.
res2 <- lmExact(x = 1:20, ny = 3, intercept = 3, slope = 2, error = 1)

## No replicates, intercept = 2 and p-value = 0.025, sigma = 3, n = 50.
## => slope = 0.063
res3 <- lmExact(x = 1:50, ny = 1, intercept = 2, pval = 0.025, error = 3)

## 5 replicates, intercept = 1, R-square = 0.85, sigma = 2, n = 10.
## => slope = 0.117
res4 <- lmExact(x = 1:10, ny = 5, intercept = 1, rsq = 0.85, error = 2)

## Heteroscedastic (magnitude-dependent) noise.
error <- sapply(1:20, function(x) rnorm(3, 0, x/10))
res5 <- lmExact(x = 1:20, ny = 3, intercept = 1, slope = 0.2,
                error = error)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

