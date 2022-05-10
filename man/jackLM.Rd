\name{jackLM}
\alias{jackLM}
\encoding{latin1}

\title{Jackknife linear model according to Quenouille 1956}

\description{
Jackknife (Leave-One-Out) method for parameter estimation and confidence interval of a linear model. 
}

\usage{
jackLM(formula, data = NULL, alpha = 0.05) 
}

\arguments{
\item{formula}{a formula of type \code{y ~ x} for the linear model.}
\item{data}{an optional data frame, list or environment containing the variables in the model.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
}

\details{
For all (\eqn{x_i, y_i}) datapoints, a linear model is created by leaving out each entry successively, \eqn{Y_{-i} = X_{-i}\beta + \varepsilon}. Pseudovalues from obtained and original coefficients are then created, \eqn{P_{-i} = (N \cdot \beta) - ((N - 1) *  \beta_{-i})}, from which the average \eqn{\overline{P_{-i}}} and standard error \eqn{\frac{\sigma}{\sqrt N}} is calculated to obtain the classical confidence interval \eqn{\overline{X}_n \pm t_{\alpha,\nu}\frac{S_n}{\sqrt{n}}}. 
}

\value{
A dataframe containg the estimated coefficients, their standard error, lower an upper confidence values and \emph{p}-values.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
Notes on bias in estimation.\cr
Quenouille MH.\cr
\emph{Biometrika}, \bold{43}, 1956, 353-36l.
}

\examples{
## Example #1 with single influencers and insignificant model (p = 0.115).
## Jackknife estimates are robust w.r.t. outlier #18.
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM1 <- lm(b ~ a)
jackLM(LM1)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

