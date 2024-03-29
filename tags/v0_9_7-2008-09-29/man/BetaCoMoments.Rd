\name{BetaCo-Moments}
\alias{BetaCoMoments}
\alias{BetaCoVariance}
\alias{BetaCoSkewness}
\alias{BetaCoKurtosis}
\alias{SystematicSkewness}
\alias{SystematicKurtosis}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to calculate systematic or beta co-moments of return series }
\description{
    calculate higher co-moment betas, which are sometimes refered to as 'systematic' skewness or kurtosis
}
\usage{
BetaCoVariance(Ra, Rb)
BetaCoSkewness(Ra, Rb)
BetaCoKurtosis(Ra, Rb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ra}{ a vector, matrix, data frame, timeSeries or zoo object of asset returns }
  \item{Rb}{ a vector, matrix, data frame, timeSeries or zoo object of index, benchmark, or secondary asset returns to compare against }
}
\details{
The co-moments do not allow the marginal impact of an asset on a portfolio to be directly measured.  Instead, Martellini and Zieman (2007) develop a framework that assesses the potential diversification of an asset relative to a portfolio. They use higher moment betas to estimate how much risk will be impacted by adding an asset to a portfolio may have, in terms of symetric risk (i.e., volatility), in asymmetry risk (i.e., skewness), and extreme risks (i.e. kurtosis). That allows them to show that adding an asset to a portfolio (or benchmark) will reduce the portfolio's variance to be reduced if the second-order beta of the asset with respect to the portfolio is less than one.  They develop the same concepts for the third and fourth order moments.  The authors offer these higher moment betas as a measure of the diversification potential of an asset.

Higher moment betas are defined as proportional to the derivative of the covariance, coskewness and cokurtosis of the second, third and fourth portfolio moment with respect to the portfolio weights. The beta co-variance is calculated as:
\begin{equation}
    \beta^{(2)}_{a,b} = \frac{CoV(R_a,R_b)}{\mu^{(2)}(R_b)}
\end{equation}
Beta co-skewness is given as:
\begin{equation}
    \beta^{(3)}_{a,b} = \frac{CoS(R_a,R_b)}{\mu^{(3)}(R_b)}
\end{equation}
Beta co-kurtosis is:
\begin{equation}
    \beta^{(4)}_{a,b} = \frac{CoK(R_a,R_b)}{\mu^{(4)}(R_b)}
\end{equation}
where the $n$-th centered moment is calculated as 
\begin{equation}
 \mu^{(n)}(R) = E\lbrack(R-E(R))^n\rbrack
\end{equation}

A beta is greater than 1 indicates that no diversification benefits should be expected from the introduction of that asset into the portfolio. Conversely, a that is less than 1 indicates that adding the new asset should reduce the resulting portfolio's volatility and kurtosis, and to an increase in skewness. More specifically, the lower the beta the higher the diversification effect on normal risk (i.e. volatility). Similarly, since extreme risks are generally characterised by negative skewness and positive kurtosis, the lower the , the higher the diversification effect on extreme risks (as reflected in Modified Value-at-Risk or ER).

The addition of a small fraction of a new asset to a portfolio leads to a decrease in the portfolio's second moment (respectively, an increase in the portfolio's third moment and a decrease in the portfolio's fourth moment) if and only if the second moment (respectively, the third moment and fourth moment) beta is less than 1 (see Martellini and Ziemann (2007) for more details).  For skewness, the interpretation is slightly more involved.  If the skewness is negative, we would expect an increase in portfolio skewness when the third moment beta is lower than 1. When the skewness of the portfolio is positive, then the condition is that the third moment beta is greater than, as opposed to lower than, one.
}
\value{
    numeric value of the normalized co-moment
}
\references{ 

Boudt, Kris, Brian G. Peterson, and Christophe Croux. 2008. Estimation and Decomposition
  of Downside Risk for Portfolios with Non-Normal Returns. Journal of Risk. Winter.

Martellini, Lionel, and Volker Ziemann. 2007. Improved Forecasts of Higher-Order Comoments and Implications for Portfolio Selection. EDHEC Risk and Asset Management Research Centre working paper.
}
\author{ Kris Boudt, Peter Carl, Brian Peterson }
\note{ 

}
\seealso{ \code{\link{CoMoments}} }
\examples{
data(managers)

BetaCoVariance(managers[, "HAM2", drop=FALSE], managers[, "SP500.TR", drop=FALSE])
BetaCoSkewness(managers[, "HAM2", drop=FALSE], managers[, "SP500.TR", drop=FALSE])
BetaCoKurtosis(managers[, "HAM2", drop=FALSE], managers[, "SP500.TR", drop=FALSE])


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ multivariate }
\keyword{ distribution }
\keyword{ models }
