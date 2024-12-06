\name{Influence plots}
\alias{inflPlot}
\alias{pvalPlot}
\encoding{latin1}

\title{Two diagnostic plots for checking p-value influencers}

\description{
Two different plot types that visualize \emph{p}-value influencers.\cr

1. \code{inflPlot}: plots the linear regression, marks the reverser(s) in darkred and displays trend lines for the full and leave-reversers-out data set (black and darkred, respectively).\cr
2. \code{pvalPlot}: plots the \emph{p}-values for each leave-one-out data point and displays the (log) \emph{p}-values as an index plot with reverser points in darkred, together with the \eqn{\alpha}-border as defined in \code{\link{lmInfl}} and the original models' \emph{p}-value.\cr
}

\usage{
inflPlot(infl, measure, ...) 
pvalPlot(infl, ...) 
}

\arguments{
\item{infl}{an object obtained from \code{\link{lmInfl}}.}
\item{measure}{which influence measure to use, see 'Details'.}
\item{...}{other plotting parameters.}
}

\details{
The influence \code{measure}'s to use are those listed in \code{\link{lmInfl}}, with the following syntax:\cr
\code{"dfb.Slope", "dffit", "cov.r", "cook.d", "hat",  "sR", "hadi", "cdr", "Si"}.
}

\value{
The corresponding plot.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
Regression diagnostics: Identifying influential data and sources of collinearity.\cr
Belsley DA, Kuh E, Welsch RE.\cr
John Wiley, New York (2004).

Applied Regression Analysis: A Research Tool.\cr
Rawlings JO, Pantula SG, Dickey DA.\cr
Springer; 2nd Corrected ed. 1998. Corr. 2nd printing 2001.

Applied Regression Analysis and Generalized Linear Models.\cr
Fox J.\cr
SAGE Publishing, 3rd ed, 2016.

Residuals and Influence in Regression.\cr
Cook RD & Weisberg S.\cr
Chapman & Hall, 1st ed, New York, USA (1982).
}

\examples{
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM1 <- lm(b ~ a)
res1 <- lmInfl(LM1) 
inflPlot(res1)
pvalPlot(res1)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

