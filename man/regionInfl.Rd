\name{regionInfl}
\alias{regionInfl}
\encoding{latin1}

\title{Identify regions of significance reversal and influence measure threshold}

\description{
Identifies regions of an (univariate) linear model in which a future data point would result in either\cr
a) significance reversal, or\cr
b) any selected influence measure as given in \code{crit} exceed its threshold value.\cr
This is intended mainly for visual/didactical purposes.
}

\usage{regionInfl(model, div.x = 20, div.y = 20, grid = TRUE, pred.int = TRUE,
       crit = c("P", "dfb.Slope", "dffit", "cov.r", "cook.d", "hat", "hadi",
       "sR", "cdr", "Si"), cex.grid = 0.5, alpha = 0.05, xlim = NULL, ylim = NULL, ...)
}

\arguments{
\item{model}{the linear model of class \code{\link{lm}}.}
\item{div.x}{the number of grid division for the \emph{x}-axis.}
\item{div.y}{the number of grid division for the \emph{y}-axis.}
\item{grid}{logical. Show the grid lines on the plot or not.}
\item{pred.int}{logical. Show the 95\% prediction interval on the plot or not.}
\item{crit}{the criterion to use. Either \code{"P"} for significance reversal or any of the influence measures given there.}
\item{cex.grid}{size of the grid points.}
\item{alpha}{the \eqn{\alpha}-level to be set as threshold.}
\item{xlim}{similar to \code{\link{xlim}}, a 2-element vector for the \emph{x}-axis limits, overrides \code{fac.x}.}
\item{ylim}{similar to \code{\link{ylim}}, a 2-element vector for the \emph{y}-axis limits, overrides \code{fac.y}.}
\item{...}{other parameters to be supplied to \code{\link{plot}} or \code{\link{lmInfl}}.}
}

\details{
For a given linear model \eqn{y_i = \beta_0 + \beta_1 x_i + \varepsilon}, each \eqn{(a, b)} pair from a grid of values \eqn{(a_1 \ldots a_j, b_1 \ldots b_k)} is added to the data, and an updated model \eqn{(y_i, b_k) = \beta_0 + \beta_1 (x_i, a_j) + \varepsilon} is created. If the updated model's \eqn{p \leq \alpha} or any of the influence measures does not exceed its published threshold, it is plotted in green, otherwise in orange. If \code{outlier = TRUE}, a possible reverser is eliminated prior to analysis but visualized in the plot.  
}

\value{
A plot with the regions marked in orange or green, and the grid matrix (\code{grid}) including the criterion outcome in 1 (green) or 0 (orange). 
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Model with p = 0.014 
set.seed(7)
N <- 20
x <- runif(N, 1, 100)
y <- 0.05 * x + rnorm(N, 0, 2)
LM1 <- lm(y ~ x)
summary(LM1)
regionInfl(LM1, crit = "P", div.x = 20, div.y = 20, cex.grid = 1, 
           xlim = c(-20, 120), ylim = c(-5, 10))
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

