\name{lmInfl}
\alias{lmInfl}
\encoding{latin1}

\title{Checks and analyzes leave-one-out (LOO) p-values and a variety of influence measures in linear regression}

\description{
This function calculates leave-one-out (LOO) \emph{p}-values for all data points and identifies those resulting in "significance reversal", i.e. in the \emph{p}-value of the model's slope traversing the user-defined \eqn{\alpha}-level. It also extends the classical influence measures from \code{\link{influence.measures}} with a few newer ones (\emph{e.g}, 'Hadi's measure', 'Coefficient of determination ratio' and 'Pena's Si') within an output format where each outlier is marked when exceeding the measure's specific threshold, as defined in the literature. Belsley, Kuh & Welsch's \emph{dfstat} criterion is also included. 
}

\usage{
lmInfl(model, alpha = 0.05, cutoff = c("BKW", "R"), verbose = TRUE, ...) 
}

\arguments{
\item{model}{the linear model of class \code{\link{lm}}.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{cutoff}{use the cutoff-values from \code{B}elsley, \code{K}uh & \code{W}elsch or the \code{R}-internal ones. See 'Details'.}
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

Cut-off values for the different influence measures are per default (\code{cutoff = "BKW"}) those defined in Belsley, Kuh & Welsch (1980) and additional literature.\cr\cr
\bold{dfbeta slope}: \eqn{| \Delta\beta1_i | > 2/\sqrt{n}} (page 28)\cr
\bold{dffits}: \eqn{| \mathrm{dffits}_i | > 2\sqrt{2/n}} (page 28)\cr
\bold{covratio}: \eqn{|\mathrm{covr}_i - 1| > 3k/n} (page 23)\cr
\bold{Cook's D}: \eqn{D_i > Q_F(0.5, k, n - k)} (Cook & Weisberg, 1982)\cr
\bold{leverage}: \eqn{h_{ii} > 2k/n} (page 17)\cr
\bold{studentized residual}: \eqn{t_i > Q_t(0.975, n - k - 1)} (page 20)\cr

If (\code{cutoff = "R"}), the criteria from \code{\link{influence.measures}} are employed:\cr\cr
\bold{dfbeta slope}: \eqn{| \Delta\beta1_i | > 1}\cr
\bold{dffits}: \eqn{| \mathrm{dffits}_i | > 3\sqrt{(k/(n - k))}}\cr
\bold{covratio}: \eqn{|1 - \mathrm{covr}_i| > 3k/(n - k)}\cr
\bold{Cook's D}: \eqn{D_i > Q_F(0.5, k, n - k)}\cr
\bold{leverage}: \eqn{h_{ii} > 3k/n}\cr 

The influence output also includes the following more "recent" measures:\cr
\bold{Hadi's measure} (column "hadi"):
\deqn{H_i^2 = \frac{h_{ii}}{1 - h_{ii}} + \frac{p}{1 - h_{ii}}\frac{d_i^2}{(1-d_i^2)}}
where \eqn{h_{ii}} are the diagonals of the hat matrix (leverages), \eqn{p = 2} in univariate linear regression and \eqn{d_i = e_i/\sqrt{\rm{SSE}}}, and threshold value \eqn{\mathrm{Med}(H_i^2) + 2 \cdot \mathrm{MAD}(H_i^2)}.

\bold{Coefficient of Determination Ratio} (column "cdr"):
\deqn{\mathrm{CDR}_i = \frac{R_{(i)}^2}{R^2}}
with \eqn{R_{(i)}^2} being the coefficient of determination without value \emph{i}, and threshold 
\deqn{\frac{B_{\alpha,p/2,(n-p-2)/2}}{B_{\alpha,p/2,(n-p-1)/2}}}

\bold{Pena's Si} (column "Si"):
\deqn{S_i = \frac{\mathbf{s}'_i\mathbf{s}_i}{p\widehat{\mathrm{var}}(\hat{y}_i)}}
where \eqn{\mathbf{s_i}} is the vector of each fitted value from the original model, \eqn{\hat{y}_i}, subtracted with all fitted values after 1-deletion, \eqn{\hat{y}_i - \hat{y}_{i(-1)}, \cdots, \hat{y}_i - \hat{y}_{i(-n)}}, \eqn{p} = number of parameters, and \eqn{\widehat{\mathrm{var}}(\hat{y}_i) = s^2h_{ii}}, \eqn{s^2 = (\mathbf{e}'\mathbf{e})/(n - p)}, \eqn{\mathbf{e}} being the residuals. In this package, a cutoff value of 0.9 is used, as the published criterion of \eqn{|\mathbf{S_i} - \mathrm{Med}(\mathbf{S})| \ge 4.5\mathrm{MAD}(\mathbf{S})} seemed too conservative. Results from this function were verified by Prof. Daniel Pena through personal communication.
}

\value{
A list with the following items:\cr
\item{origModel}{the original model with all data points.} 
\item{finalModels}{a list of final models with the influencer(s) removed.}
\item{infl}{a matrix with the original data, classical \code{\link{influence.measures}}, studentized residuals, leverages, dfstat, LOO-\emph{p}-values, LOO-slopes/intercepts and their \eqn{\Delta}'s, LOO-standard errors and \eqn{R^2}s. Influence measures that exceed their specific threshold - see \code{\link{inflPlot}} - will be marked with asterisks.}
\item{raw}{same as \code{infl}, but with pure numeric data.}
\item{sel}{a vector with the influencers' indices.}
\item{alpha}{the selected \eqn{\alpha}-level.}
\item{origP}{the original model's \emph{p}-value.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{For dfstat / dfstud :}\cr
Regression diagnostics: Identifying influential data and sources of collinearity.\cr
Belsley DA, Kuh E, Welsch RE.\cr
John Wiley, New York, USA (2004).

Econometrics, 5ed.\cr
Baltagi B.\cr
Springer-Verlag Berlin, Germany (2011).

Growth regressions and what the textbooks don't tell you.\cr
Temple J.\cr
\emph{Bull Econom Res}, \bold{52}, 2000, 181-205.

Robust Regression and Outlier Detection.\cr
Rousseeuw PJ & Leroy AM.\cr
John Wiley & Sons, New York, NY (1987).\cr

\bold{Hadi's measure:}\cr
A new measure of overall potential influence in linear regression.\cr
Hadi AS.\cr
\emph{Comp Stat & Data Anal}, \bold{14}, 1992, 1-27.\cr

\bold{Coefficient of determination ratio:}\cr
On the detection of influential outliers in linear regression analysis.\cr
Zakaria A, Howard NK, Nkansah BK.\cr
\emph{Am J Theor Appl Stat}, \bold{3}, 2014, 100-106.

On the Coefficient of Determination Ratio for Detecting Influential Outliers in Linear Regression Analysis.\cr
Zakaria A, Gordor BK, Nkansah BK.\cr
\emph{Am J Theor Appl Stat}, \bold{11}, 2022, 27-35.\cr

\bold{Pena's measure:}\cr
A New Statistic for Influence in Linear Regression.\cr
Pena D.\cr
\emph{Technometrics}, \bold{47}, 2005, 1-12.\cr
}

\examples{
## Example #1 with single influencer and significant model (p = 0.0089).
## Removal of #21 results in p = 0.115!
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
a <- c(a, 25); b <- c(b, 10)
LM1 <- lm(b ~ a)
lmInfl(LM1) 

## Example #2 with single influencer and insignificant model (p = 0.115).
## Removal of #18 results in p = 0.0227!
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM2 <- lm(b ~ a)
lmInfl(LM2) 

## Example #3 with multiple influencers and significant model (p = 0.0269).
## Removal of #2, #17, #18 or #20 results in crossing p = 0.05!
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM3 <- lm(b ~ a)
lmInfl(LM3) 

## Large Example #4 with top 10 influencers and significant model (p = 6.72E-8).
## Not possible to achieve a crossing of alpha with any point despite strong noise.
set.seed(123)
a <- 1:100
b <- 5 + 0.08 * a + rnorm(100, 0, 5)
LM4 <- lm(b ~ a)
lmInfl(LM4) 
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

