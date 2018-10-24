\name{lmMult}
\alias{lmMult}
\encoding{latin1}

\title{Checks and analyzes leave-multiple-out (LMO) p-values in linear regression}

\description{
This function calculates leave-multiple-out (LMO) \emph{p}-values for an increasing number of datapoints and identifies those resulting in "significance reversal" of the model, i.e. in the slope's \emph{p}-value traversing the user-defined \eqn{\alpha}-level.
}

\usage{
lmMult(model, max = 5, n = 10000, alpha = 0.05, 
       method = c("pearson", "spearman"), verbose = TRUE) 
}

\arguments{
\item{model}{the linear model of class \code{\link{lm}}.}
\item{max}{the maximum number of points to eliminate.}
\item{n}{the number of samples to draw for each 1...\code{max}.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{method}{select either parametric (\code{"pearson"}) or rank-based (\code{"spearman"}) statistics.}
\item{verbose}{logical. If \code{TRUE}, results for each 1...\code{max} will be printed to the console.}
}

\details{
The algorithm\cr
1) calculates the \emph{p}-value of the full model (all data points),\cr
2) calculates a LMO-\emph{p}-value for all \code{n} sampled groups of 1...\code{max} points removed,\cr
3) checks for significance reversal in the resulting model and\cr
4) returns all \code{n} samples and the corresponding \emph{p}-values.\cr
}

\value{
A list with the following items:\cr
\item{sample}{a matrix with all \code{max * n} iterations, where a 1 indicates the left-out sample(s), as well as the corrsponding \emph{p}-values and group.}
\item{stat}{for each 1...\code{max} LMO's, the percentage of model significance reversers.}
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Example with single influencers and insignificant model (p = 0.115).
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM1 <- lm(b ~ a)
res1 <- lmMult(LM1)
multPlot(res1)
stability(res1)

## Large example with 100 data points and highly significant model (p = 6.72E-8).
## No significance reversal up to the elimination of 20 points.
\dontrun{
set.seed(123)
a <- 1:100
b <- 5 + 0.08 * a + rnorm(100, 0, 5)
LM2 <- lm(b ~ a)
res2 <- lmMult(LM2, max = 20)
multPlot(res2)
stability(res2)
}
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

