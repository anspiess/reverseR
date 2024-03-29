\name{Influence plots}
\alias{lmPlot}
\alias{pvalPlot}
\alias{inflPlot}
\alias{slsePlot}
\alias{threshPlot}
\alias{stabPlot}
\encoding{latin1}

\title{Several diagnostic plots for checking p-value influencers}

\description{
Seven different plot types that visualize \emph{p}-value influencers.\cr

1. \code{lmPlot}: plots the linear regression, marks the influencer(s) in red and displays trend lines for the full and leave-one-out (LOO) data set (black and red, respectively).\cr
2. \code{pvalPlot}: plots the \emph{p}-values for each LOO data point and displays the values as a full model/LOO model plot, together with the \code{alpha} border as defined in \code{\link{lmInfl}}.\cr
3. \code{inflPlot}: plots \code{\link{dfbeta}} for slope, \code{\link{dffits}}, \code{\link{covratio}}, \code{\link{cooks.distance}}, leverage (\code{\link{hatvalues}}), studentized residuals (\code{\link{rstudent}}) and Hadi's measure against the \eqn{\Delta}\emph{p}-value. Herewith, changes in these six parameters can be compared to the effect on the corresponding drop/rise in \emph{p}-value. The plots include vertical boundaries for threshold values as defined in the literature under 'References'.\cr
4. \code{slsePlot}: plots all LOO-slopes and their standard errors together with the corresponding original model values and a \emph{t}-value border as calculated by \eqn{\mathit{Q_t}(1 - \frac{\alpha}{2}, n-2)}. LOO of points on the right of this border result in a significant model, and \emph{vice versa}.\cr
5. \code{threshPlot}: plots the output of \code{\link{lmThresh}}, i.e. the regression plot including confidence/prediction intervals, as well as for each response value \eqn{y_i} the region in which the model is significant (green). This is tested for either i) \eqn{y_i} that are shifted into this region (\code{newobs = FALSE} in \code{\link{lmThresh}}) or ii) when a new observation \eqn{y2_i} is added (\code{newobs = TRUE} in \code{\link{lmThresh}}). In the latter case, it is informative if this region resides within the prediction interval (dashed line), indicating that a future additional measurement at \eqn{x_i} might reverse the significance statement.\cr
6. \code{stabPlot}: for single (to be selected) response values from the output of \code{\link{lmThresh}}, this function displays the region of significance reversal within the surrounding prediction interval. The probability of a either shifting the response value (if \code{lmThresh(..., newobs = FALSE)}) or of including a future (measurement) point (if \code{lmThresh(..., newobs = TRUE)}) to reverse the significance is shown as the integral between the "end of significance region" (eosr) and the nearest prediction interval boundary.

\bold{NOTE}: The visual display should always be supplemented with the corresponding \code{\link{stability}} analysis.
}

\usage{
lmPlot(infl, ...) 
pvalPlot(infl, ...) 
inflPlot(infl, ...)
slsePlot(infl, ...)
threshPlot(thresh, bands = FALSE, ...)
stabPlot(stab, which = NULL, ...)
}

\arguments{
\item{infl}{an object obtained from \code{\link{lmInfl}}.}
\item{thresh}{an object obtained from \code{\link{lmThresh}}.}
\item{stab}{an object obtained from using \code{\link{stability}} on an \code{\link{lmThresh}} output.}
\item{bands}{logical. If \code{TRUE}, plots the confidence and prediction bands.}
\item{which}{which response value should be shown in \code{stabPlot}?}
\item{...}{other plotting parameters.}
}

\value{
The corresponding plot.
}

\note{
Cut-off values for the different influence measures are those defined in Belsley, Kuh E & Welsch (1980):\cr\cr
\bold{dfbeta slope}: \eqn{| \Delta\beta1_i | > 2/\sqrt{n}} (page 28)\cr
\bold{dffits}: \eqn{| \mathrm{dffits}_i | > 2\sqrt{2/n}} (page 28)\cr
\bold{covratio}: \eqn{|\mathrm{covr}_i - 1| > 3k/n} (page 23)\cr
\bold{Cook's D}: \eqn{D_i > Q_F(0.5, k, n - k)} (Cook & Weisberg, 1982)\cr
\bold{leverage}: \eqn{h_{ii} > 2k/n} (page 17)\cr
\bold{studentized residual}: \eqn{t_i > Q_t(0.975, n - k - 1)} (page 20)\cr
\bold{Hadi's measure}: \eqn{H_i^2 > \mathrm{Med}(H_i^2) + 2 \cdot \mathrm{MAD}(H_i^2)} (Hadi 1992)
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

A new measure of overall potential influence in linear regression.\cr
Hadi AS.\cr
\emph{Comp Stat & Data Anal}, \bold{14}, 1992, 1-27.
}

\examples{
## See Examples in 'lmInfl' and 'lmThresh'.
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

