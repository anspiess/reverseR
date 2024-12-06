\name{pcomp}
\alias{pcomp}
\encoding{latin1}

\title{Calculates linear regression p-values from a variety of robust regression methods}

\description{
This function calculates \emph{p}-values from a variety of methods, specifically:\cr
1) standard linear model\cr
2) standard linear model with highest \emph{p}-influencer removed\cr
3) robust regression with MM-estimators\cr
4) Theil-Sen regression\cr
5) least absolute deviations regression\cr
6) quantile regression\cr
7) weighted regression with isolation forest scores as inverse weights\cr
8) bootstrap linear model, see \code{\link{bootLM}}\cr
9) jackknife linear model, see \code{\link{jackLM}}
}

\usage{
pcomp(x, y = NULL, R = 1000, alpha = 0.05, ...)
}

\arguments{
\item{x}{either a linear model of class \code{\link{lm}} or the regressions \emph{x}-values.}
\item{y}{the optional \emph{y}-values.}
\item{R}{the number of bootstrap resamples, see \code{\link{bootLM}}.}
\item{alpha}{the \eqn{\alpha}-level for \code{\link{lmInfl}}.}
\item{...}{further arguments to be passed to downstream methods.}
}

\details{
This function is meant to provide a swift overview on the sensitivity of the \emph{p}-values to different (mostly robust) linear regression methods, which correlates to a large extent with the presence of influential / outlying data points, see 'Examples'.
}

\value{
A vector of \emph{p}-values from the above mentioned ten methods, in that order.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
Robust Regression and Outlier Detection.\cr
Rousseeuw PJ & Leroy AM.\cr
1ed (1987), Wiley (NJ, USA).

A rank-invariant method of linear and polynomial regression analysis.\cr
Theil H.\cr 
\emph{I. Nederl. Akad. Wetensch. Proc}, \bold{53}, 1950, 386-392.

Estimates of the regression coefficient based on Kendall's tau.\cr
Sen PK.\cr
\emph{J Am Stat Assoc}, \bold{63}, 1968, 1379-1389.

Least absolute deviations estimation via the EM algorithm.\cr
Phillips RF.\cr
\emph{Statistics and Computing}, \bold{12}, 2002, 281-285. 

Quantile Regression.\cr
Koenker R.\cr
Cambridge University Press, Cambridge, New York (2005).

Isolation-based anomaly detection.\cr
Liu FT, Ting KM, Zhou ZH.\cr
\emph{ACM Transactions on Knowledge Discovery from Data}, \bold{6.1}, 2012, 3.
}

\examples{
## Example with influencer
## => a few methods indicate significant 
## downward drop of the p-value
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
pcomp(a, b) 
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

