\name{simInfl}
\alias{simInfl}
\encoding{latin1}

\title{Simulates significance reversals and calculates their influence parameters}

\description{
This function simulates linear regressions and stores the parameters and influence measures of all simulations that resulted in LOO significance reversal, developed for research purposes.
}

\usage{
simInfl(x = 1:10, slope = 0.02, intercept = 1, error = 0.05, nrev = 1000, ...) 
}

\arguments{
\item{x}{the \eqn{x} values to be supplied to \code{\link{lmExact}}.}
\item{slope}{the slope \eqn{\beta_1} to be supplied to \code{\link{lmExact}}.}
\item{intercept}{the intercept \eqn{\beta_0} to be supplied to \code{\link{lmExact}}.}
\item{error}{the \eqn{\varepsilon} value to be supplied to \code{\link{lmExact}}.}
\item{nrev}{the number of desired significance reversals.}
\item{...}{other parameters to \code{\link{lmExact}} and \code{\link{lmInfl}}.}
}

\details{
Loops over an undefined number of EXACT regressions (\code{\link{lmExact}}) with incrementing random seeds, stores all models and in case of significance reversal, parameters and influence measures (\code{\link{lmInfl}}). The simulation terminates when \code{nrev} reversals are counted.
}

\value{
A list with the following two items:
\item{models}{the linear models of all reversals.}
\item{mat}{the stored matrix with the resulting parameters and influence measures for all \code{nrev} reversals.}
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
\dontrun{
## Example with slight slope, intercept = 0.5 and 100 reversals.
res <- simInfl(x = 1:10, intercept = 0.5, slope = 0.02, error = 0.05, nrev = 100)

## Plot Cook's D versus delta-P values
## and insert common cut-off.
plot(res$mat[, "cook.d"], res$mat[, "dP"], pch = 16, cex = 0.5,
     xlab = "Cook's D", ylab = "delta-P")
thresh <- qf(0.5, 2, 8)  # threshold value for Qf(0.5, npar, df)
abline(v = thresh, col = "darkred", lwd = 2)  

## Plot dfbeta slope versus delta-P values
## and insert common cut-off.
plot(res$mat[, "dfb.Slope"], res$mat[, "dP"], pch = 16, cex = 0.5,
     xlab = "dfbeta Slope", ylab = "delta-P")
thresh <- 2/sqrt(10)  # 2/sqrt(N)
abline(v = thresh, col = "darkred", lwd = 2)  

## Plot dffits versus delta-P values
## and insert common cut-off.
plot(abs(res$mat[, "dffit"]), res$mat[, "dP"], pch = 16, cex = 0.5,
     xlab = "dffits", ylab = "delta-P")
thresh <- 2 * sqrt(2/10)  # 2 * sqrt(nPar/N)
abline(v = thresh, col = "darkred", lwd = 2)  
}
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

