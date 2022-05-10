\name{bootLM}
\alias{bootLM}
\encoding{latin1}

\title{Nonparametric bootstrap linear model}

\description{
Nonparametric bootstrap (sampling cases with replacement) method for parameter estimation and confidence interval of a linear model. 
}

\usage{
bootLM(formula, data = NULL, R = 10000, alpha = 0.05) 
}

\arguments{
\item{formula}{a formula of type \code{y ~ x} for the linear model.}
\item{data}{an optional data frame, list or environment containing the variables in the model.}
\item{R}{number of bootstrap samples.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
}

\details{
For all (\eqn{x_i, y_i}) datapoints, linear models are created by sampling \code{R} times - with replacement - from \eqn{n \in \{1 \ldots N\}} and building models \eqn{Y_n = X_n\beta + \varepsilon}. This is also known as the .632-bootstrap, because the samples will, on average, contain \eqn{1 - e^{-1} = 0.632} unique elements. 
Parameter estimates are obtained from each sampling, from which the average \eqn{\overline{P_{n}}} and standard error \eqn{\frac{\sigma}{\sqrt n}} is calculated as well as a quantile based confidence interval. \emph{p}-values are calculated through inversion of the confidence interval (\code{\link{boot.pval}}). 
}

\value{
A dataframe containing the estimated coefficients, their standard error, lower an upper confidence values and \emph{p}-values.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
An Introduction to the Bootstrap.\cr
Efron B, Tibshirani R.\cr
Chapman & Hall (1993).

The Bootstrap and Edgeworth Expansion.\cr
Hall P.\cr 
Springer, New York (1992).

Modern Statistics with R.\cr
Thulin M.\cr
Eos Chasma Press, Uppsala (2021).
}

\examples{
## Example #1 with single influencers and insignificant model (p = 0.115).
## Jackknife estimates are robust w.r.t. outlier #18.
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM1 <- lm(b ~ a)
bootLM(LM1, R = 1000)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

