\name{lmInfl}
\alias{lmInfl}
\encoding{latin1}

\title{Checks and analyzes leave-one-out (LOO) p-values in linear regression}

\description{
This function calculates leave-one-out (LOO) \emph{p}-values for all data points and identifies those resulting in "significance reversal", i.e. in the \emph{p}-value of the model's slope traversing the user-defined \eqn{\alpha}-level.
}

\usage{
lmInfl(model, alpha = 0.05, verbose = TRUE, ...) 
}

\arguments{
\item{model}{the linear model of class \code{\link{lm}}.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{verbose}{logical. If \code{TRUE}, results are displayed on the console.}
\item{...}{other arguments to \code{\link{lm}}.}
}

\details{
The algorithm\cr
1) calculates the \emph{p}-value of the full model (all points),\cr
2) calculates a LOO-\emph{p}-value for each point removed,\cr
3) checks for significance reversal in all data points and\cr
4) returns all models as well as classical \code{\link{influence.measures}} with LOO-\emph{p}-values, \eqn{\Delta}\emph{p}-values, slopes and standard errors attached.\cr

The idea of \emph{p}-value influencers was first introduced by Belsley, Kuh & Welsch, and described as an influence measure pertaining directly to the change in \emph{t}-statistics, that will "show whether the conclusions of hypothesis testing would be affected", termed \bold{dfstat} in [1, 2, 3] or \bold{dfstud} in [4]:
\deqn{\rm{dfstat}_{ij} \equiv \frac{\hat{\beta}_j}{s\sqrt{(X'X)^{-1}_{jj}}}-\frac{\hat{\beta}_{j(i)}}{s_{(i)}\sqrt{(X'_{(i)}X_{(i)})^{-1}_{jj}}}}
where \eqn{\hat{\beta}_j} is the \emph{j}-th estimate, \emph{s} is the residual standard error, \emph{X} is the design matrix and (\emph{i}) denotes the \emph{i}-th observation deleted.\cr
\bold{dfstat}, which for the regression's slope \eqn{\beta_1} is the difference of \emph{t}-statistics 
\deqn{\Delta t = t_{\beta1} - t_{\beta1(i)} = \frac{\beta_1}{\rm{s.e.(\beta_1)}} - \frac{\beta_1(i)}{\rm{s.e.(\beta_1(i)})}}
is inextricably linked to the changes in \emph{p}-value \eqn{\Delta p}, calculated from
\deqn{\Delta p = p_{\beta1} - p_{\beta1(i)} = 2\left(1-P_t(t_{\beta1}, \nu)\right) - 2\left(1-P_t(t_{\beta1(i)} , \nu-1)\right)}
where \eqn{P_t} is the Student's \emph{t} cumulative distribution function with \eqn{\nu} degrees of freedom, and where significance reversal is attained when \eqn{\alpha \in [p_{\beta1}, p_{\beta1(i)}]}.
Interestingly, the seemingly mandatory check of the influence of single data points on statistical inference is living in oblivion: apart from [1-4], there is, to the best of our knowledge, no reference to \bold{dfstat} or \eqn{\Delta p} in current literature on influence measures.

The influence output also includes the more recent Hadi's measure (column "hadi"):
\deqn{H_i^2 = \frac{p_{ii}}{1 - p_{ii}} + \frac{k}{1 - p_{ii}}\frac{d_i^2}{(1-d_i^2)}}
where \eqn{p_{ii}} are the diagonals of the hat matrix (leverages), \eqn{k = 2} in univariate linear regression and \eqn{d_i = e_i/\sqrt{\rm{SSE}}}.
}

\value{
A list with the following items:\cr
\item{origModel}{the original model with all data points.} 
\item{finalModels}{a list of final models with the influencer(s) removed.}
\item{infl}{a matrix with the original data, classical \code{\link{influence.measures}}, studentized residuals, leverages, LOO-\emph{p}-values, LOO-slopes/intercepts and their \eqn{\Delta}'s, LOO-standard errors and \eqn{R^2}s. Influence measures that exceed their specific threshold - see \code{\link{inflPlot}} - will be marked with asterisks.}
\item{raw}{same as \code{infl}, but with pure numeric data.}
\item{sel}{a vector with the influencers' indices.}
\item{alpha}{the selected \eqn{\alpha}-level.}
\item{origP}{the original model's \emph{p}-value.}
\item{stab}{the stability measure, see \code{\link{stability}}.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{For dfstat / dfstud :}\cr
1. Regression diagnostics: Identifying influential data and sources of collinearity.\cr
Belsley DA, Kuh E, Welsch RE.\cr
John Wiley, New York, USA (2004).\cr

2. Econometrics, 5ed.\cr
Baltagi B.\cr
Springer-Verlag Berlin, Germany (2011).\cr

3. Growth regressions and what the textbooks don't tell you.\cr
Temple J.\cr
\emph{Bull Econom Res}, \bold{52}, 2000, 181-205.\cr

4. Robust Regression and Outlier Detection.\cr
Rousseeuw PJ & Leroy AM.\cr
John Wiley & Sons, New York, NY (1987).\cr

\bold{Hadi's measure:}\cr
A new measure of overall potential influence in linear regression.\cr
Hadi AS.\cr
\emph{Comp Stat & Data Anal}, \bold{14}, 1992, 1-27.
}

\examples{
## Example #1 with single influencers and insignificant model (p = 0.115).
## Removal of #18 results in p = 0.0227!
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM1 <- lm(b ~ a)
res1 <- lmInfl(LM1) 
lmPlot(res1)
pvalPlot(res1)
inflPlot(res1)
slsePlot(res1)
stability(res1)

## Example #2 with multiple influencers and significant model (p = 0.0269).
## Removal of #2, #17, #18 or #20 result in crossing p = 0.05!
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM2 <- lm(b ~ a)
res2 <- lmInfl(LM2) 
lmPlot(res2)
pvalPlot(res2)
inflPlot(res2)
slsePlot(res2)
stability(res2)

## Large Example #3 with top 10 influencers and significant model (p = 6.72E-8).
## Not possible to achieve a crossing of alpha with any point despite strong noise.
set.seed(123)
a <- 1:100
b <- 5 + 0.08 * a + rnorm(100, 0, 5)
LM3 <- lm(b ~ a)
res3 <- lmInfl(LM3) 
lmPlot(res3)
stability(res3)

## Example #4 with replicates and single influencer (p = 0.114).
## Removal of #58 results in p = 0.039.
set.seed(123)
a <- rep(1:20, each = 3)
b <- 5 + 0.08 * a + rnorm(20, 0, 2)
LM4 <- lm(b ~ a)
res4 <- lmInfl(LM4) 
lmPlot(res4)
pvalPlot(res4)
inflPlot(res4)
slsePlot(res4)
stability(res4)

## As Example #1, but with weights.
## Removal of #18 results in p = 0.04747.
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM5 <- lm(b ~ a, weights = 1:20)
res5 <- lmInfl(LM5) 
lmPlot(res5)
stability(res5)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

