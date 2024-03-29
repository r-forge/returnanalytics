\name{chart.QQPlot}
\alias{chart.QQPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a QQ chart }
\description{
Plot the return data against any theoretical distribution.
}
\usage{
chart.QQPlot(R, distribution = "norm", ylab = NULL,
                 xlab = paste(distribution, "Quantiles"), main = NULL,
                 las = par("las"), envelope = FALSE, labels = FALSE,
                 col = c(1, 4), lwd = 2, pch = 1, cex = 1,
                 line = c("quartiles", "robust", "none"),
                 elementcolor = "lightgray", ...)
}
\details{
A Quantile-Quantile (QQ) plot is a scatter plot designed to compare the data to the theoretical distributions to visually determine if the observations are likely to have come from a known population. The empirical quantiles are plotted to the y-axis, and the x-axis contains the values of the theorical model.  A 45-degree reference line is also plotted. If the empirical data come from the population with the choosen distribution, the points should fall approximately along this reference line. The larger the departure from the reference line, the greater the evidence that the data set have come from a population with a different distribution.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{ a vector, matrix, data frame, timeSeries or zoo object of asset returns }
  \item{distribution}{ root name of comparison distribution - e.g., 'norm' for
          the normal distribution; 't' for the t-distribution. See examples for other ideas. }
  \item{xlab}{ set the x-axis label, as in \code{\link{plot}} }
  \item{ylab}{ set the y-axis label, as in \code{\link{plot}} }
  \item{main}{ set the chart title, same as in \code{plot} }
  \item{las}{ set the direction of axis labels, same as in \code{plot} }
  \item{envelope}{ confidence level for point-wise confidence envelope, or
          'FALSE' for no envelope. }
  \item{labels}{ vector of point labels for interactive point identification,
          or 'FALSE' for no labels. }
  \item{col}{ color for points and lines; the default is the \emph{second} entry
          in the current color palette (see 'palette' and 'par'). }
  \item{lwd}{ set the line width, as in \code{\link{plot}} }
  \item{pch}{ symbols to use, see also \code{\link{plot}} }
  \item{cex}{ symbols to use, see also \code{\link{plot}} }
  \item{line}{ '"quartiles"' to pass a line through the quartile-pairs, or
          '"robust"' for a robust-regression line; the latter uses the
          'rlm' function in the 'MASS' package. Specifying 'line =
          "none"' suppresses the line. }
  \item{elementcolor}{ provides the color for drawing less-important chart elements, such as the box lines, axis lines, etc. replaces \code{darken} }
  \item{\dots}{ any other passthru parameters to the distribution function }
}
\value{
QQ Plot
}
\references{
  main code forked/borrowed/ported from the excellent: \\
  Fox, John (2007) \emph{car: Companion to Applied Regression} \\
  \url{http://www.r-project.org}, \url{http://socserv.socsci.mcmaster.ca/jfox/}
 }
\author{ John Fox, ported by Peter Carl }
\seealso{
  \code{\link[stats]{qqplot}} \\
  \code{\link[car]{qq.plot}}  \\
  \code{\link{plot}}
}
\examples{
library(MASS)
data(managers)
x = checkData(managers[,2, drop = FALSE], na.rm = TRUE, method = "vector")
layout(rbind(c(1,2),c(3,4)))
# Panel 1, Normal distribution
chart.QQPlot(x, main = "Normal Distribution", distribution = 'norm', envelope=0.95)
# Panel 2, Log-Normal distribution
fit = fitdistr(1+x, 'lognormal')
chart.QQPlot(1+x, main = "Log-Normal Distribution", envelope=0.95, distribution='lnorm', meanlog = fit$estimate[[1]], sdlog = fit$estimate[[2]])
# Panel 3, Skew-T distribution
library(sn)
fit = st.mle(y=x)
chart.QQPlot(x, main = "Skew T Distribution", envelope=0.95, distribution = 'st', location = fit$dp[[1]], scale = fit$dp[[2]], shape = fit$dp[[3]], df=fit$dp[[4]], log = FALSE)
#Panel 4: Stable Parietian
library(fBasics)
fit.stable = stableFit(x,doplot=FALSE)
chart.QQPlot(x, main = "Stable Paretian Distribution", envelope=0.95, distribution = 'stable', alpha = fit.stable@fit$estimate[[1]], beta = fit.stable@fit$estimate[[2]], gamma = fit.stable@fit$estimate[[3]], delta = fit.stable@fit$estimate[[4]], pm = 0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ multivariate }
\keyword{ distribution }
\keyword{ models }
\keyword{ hplot }
