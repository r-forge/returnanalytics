\name{maxDrawdown}
\alias{maxDrawdown}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ caclulate the maximum drawdown from peak equity }
\description{
To find the maximum drawdown in a return series, we need to first calculate the cumulative returns and the maximum cumulative return to that point.  Any time the cumulative returns dips below the maximum cumulative returns, it's a drawdown.  Drawdowns are measured as a percentage of that maximum cumulative return, in effect, measured from peak equity.
}
\usage{
maxDrawdown(R, geometric = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{ an xts, vector, matrix, data frame, timeSeries or zoo object of asset returns }
  \item{geometric}{ generate geometric (TRUE) or simple (FALSE) returns, default TRUE }
  \item{\dots}{ any other passthru parameters }
}
\references{
Bacon, C. \emph{Practical Portfolio Performance Measurement and Attribution}. Wiley. 2004. p. 88 \cr

}
\author{ Peter Carl }
\seealso{
    \code{\link{findDrawdowns}} \cr
    \code{\link{sortDrawdowns}} \cr
    \code{\link{table.Drawdowns}} \cr
    \code{\link{table.DownsideRisk}} \cr
    \code{\link{chart.Drawdown}} \cr
}
\examples{
data(edhec)
maxDrawdown(edhec[,"Funds of Funds"])
data(managers)
maxDrawdown(managers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ multivariate }
\keyword{ distribution }
\keyword{ models }
