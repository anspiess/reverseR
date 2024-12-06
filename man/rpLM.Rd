\name{rpLM}
\alias{rpLM}
\encoding{latin1}

\title{Calculates the 'replication probability of significance' of an 'lm' object}

\description{
This function uses a bootstrap approach to calculate the \emph{replication probability} of significance, which answers the question "if we repeat this linear regression under identical conditions (similar sample size, similar residual variance), what is the probability of observing significance (or non-significance) similar to the original data?". 
}

\usage{
rpLM(model, alpha = 0.05, R = 10000, plot = TRUE, verbose = TRUE, ...) 
}

\arguments{
\item{model}{a linear model of class \code{\link{lm}}.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{R}{the number of bootstrap resamples, see \code{\link{bootLM}}.}
\item{plot}{logical. If \code{TRUE}, a stripchart of the bootstrap \emph{P}-values, the original \emph{P}-value and the \eqn{\alpha}-level is displayed.}
\item{verbose}{logical. If \code{TRUE}, the analysis steps are written to the console.}
\item{...}{other parameters to be supplied to \code{\link{bootLM}}.}
}

\details{
The approach here is along the lines of Boos & Stefanski (2011), which investigated the replication probability of the \emph{P}-value, as opposed to the works of 
Goodman (1992), Shao & Chow (2002) and Miller (2009), where the effect size is used. In our context, for a given linear model and using a bootstrap approach, the \emph{replication probability} is the proportion of bootstrap \emph{P}-values with \eqn{\tilde{P} \leq \alpha} when the original model is significant, or  \eqn{\tilde{P} > \alpha} when not. Hence, we employ the bootstrap to assess the sampling variability of the \emph{P}-value, not the sampling variability of the \emph{P}-value under \eqn{H_0}, as is common, thereby preserving the non-null property of the data.

Bootstrap results are obtained from non-parametric cases bootstrapping ("np.cases"), non-parametric residuals bootstrapping ("np.resid") and parametric residuals bootstrapping ("p.resid"), see \code{\link{bootLM}}.
}

\value{
A vector with the three different bootstrap results as described above. 
}

\author{
Andrej-Nikolai Spiess
}   

\references{
Ecological Models and Data in R.\cr
Chapter 5: Stochastic simulation and power analysis.\cr
Benjamin M. Bolker.\cr
Princeton University Press (2008).

P-Value Precision and Reproducibility.\cr
Boos DD & Stefanski LA.\cr
\emph{Am Stat}, \bold{65}, 2011, 213-212.

A comment on replication, p-values and evidence.\cr
Goodman SN.\cr
\emph{Stat Med}, \bold{11}, 1992, 875-879.

Reproducibility probability in clinical trials.\cr
Shao J & Chow SC.\cr 
\emph{Stat Med}, \bold{21}, 2002, 1727-1742.

What is the probability of replicating a statistically significant effect?\cr
Miller J.\cr
\emph{Psych Bull & Review}, \bold{16}, 2009, 617-640.
}

\examples{
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(length(a), 0, 1)
LM1 <- lm(b ~ a)
summary(LM1)
rpLM(LM1, R = 100)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

