\name{bootLM}
\alias{bootLM}
\encoding{latin1}

\title{Nonparametric/Parametric bootstrap linear model}

\description{
Nonparametric and parametric bootstrap (sampling cases, residuals or distributions with replacement) method for parameter estimation and confidence interval of a linear model. 
}

\usage{
bootLM(model, type = c("cases", "residuals", "residuals2", "parametric"), 
       R = 10000, alpha = 0.05, ret.models = FALSE) 
}

\arguments{
\item{model}{an \code{\link{lm}} model.}
\item{type}{what to bootstrap. See "Details".}
\item{R}{number of bootstrap samples.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{ret.models}{logical. If \code{TRUE}, the \code{R} models are returned as a list.}
}

\details{
If \code{type = "cases"}, for all (\eqn{x_i, y_i}) datapoints, linear models are created by sampling \code{R} times - with replacement - from \eqn{n \in \{1 \ldots N\}} and building models \eqn{Y_n = X_n\beta + \varepsilon}. This is also known as the .632-bootstrap, because the samples will, on average, contain \eqn{1 - e^{-1} = 0.632} unique elements.
If \code{type = "residuals"}, for all residuals (\eqn{r_i = y_i - \hat{y}_i}), linear models are created by sampling \code{R} times - with replacement - from \eqn{n \in (1 \ldots N)} and building models \eqn{\hat{Y}_i + r_n = X_i\beta + \varepsilon}. If \code{type = "residuals2"} is selected, scaled and centered residuals \eqn{r_n = \frac{r_i}{\sqrt{1 - h_{ii}}} - \bar{r}} according to Davison & Hinkley are used. In the \code{"parametric"} bootstrap, \eqn{n} values drawn from a normal distribution \eqn{j_n \in \mathcal{N}(0, \sigma)}, where \eqn{\sigma = \sqrt{\frac{\sum(r_i)^2}{n - p}}}, are added to the fitted values, and linear models are created \eqn{\hat{Y}_i + j_n = X_i\beta + \varepsilon}.
Parameter estimates are obtained from each sampling, from which the average \eqn{\overline{P_{n}}} and standard error \eqn{\hat{\sigma}} is calculated as well as a quantile based confidence interval. \emph{p}-values are calculated through inversion of the confidence interval. 
}

\value{
A dataframe containing the estimated coefficients, their standard error, lower an upper confidence values and \emph{p}-values. If \code{ret.models = TRUE} a list with all \code{R} models is returned.
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

Bootstrap methods and their application.\cr
Davison AC, Hinkley DV.\cr
Cambridge University Press (1997).
}

\examples{
## Example with single influencer (#18) and insignificant model (p = 0.115),
## using case bootstrap.
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM <- lm(b ~ a)
bootLM(LM, R = 100)

## using residuals bootstrap.
bootLM(LM, R = 100, type = "residuals")
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}