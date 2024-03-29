\name{Omega}
\alias{Omega}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calculate Omega for a return series }
\description{
Keating and Shadwick (2002) proposed Omega (referred to as Gamma in their original paper) as a way to capture all of the higher moments of the returns distribution.
}
\usage{
Omega(R, L = 0, method = c("simple", "interp", "binomial", "blackscholes"), output = c("point", "full"), rf = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{ a vector, matrix, data frame, or timeSeries of returns }
  \item{L}{ L is the loss threshold that can be specified as zero, return from a benchmark index, or an absolute rate of return - any specified level }
  \item{method}{ one of: simple, interp, binomial, blackscholes }
  \item{output}{ one of: point (in time), or full (distribution of Omega) }
  \item{rf}{ risk free rate, in same period as your returns }
}
\details{
Mathematically, Omega is:
  integral[L to b](1 - F(r))dr / integral[a to L](F(r))dr

where the cumulative distribution F is defined on the interval (a,b). L is the loss threshold that can be specified as zero, return from a benchmark index, or an absolute rate of return - any specified level. When comparing alternatives using Omega, L should be common.

Input data can be transformed prior to calculation, which may be useful for introducing risk aversion.

This function returns a vector of Omega, useful for plotting.  The steeper, the less risky.  Above it's mean, a steeply sloped Omega also implies a very limited potential for further gain.

Omega has a value of 1 at the mean of the distribution.

Omega is sub-additive.  The ratio is dimensionless.

Kazemi, Schneeweis, and Gupta (2003), in "Omega as a Performance Measure" show that Omega can be written as:
  Omega(L) = C(L)/P(L)
where C(L) is essentially the price of a European call option written on the investment and P(L) is essentially the price of a European put option written on the investment.  The maturity for both options is one period (e.g., one month) and L is the strike price of both options.

The numerator and the denominator can be expressed as:
  exp(-rf) * E[max(x - L, 0)]
  exp(-rf) * E[max(L - x, 0)]
with exp(-rf) calculating the present values of the two, where rf is the per-period riskless rate.

The first three methods implemented here focus on that observation. The first method takes the simplification described above.  The second uses the Black-Scholes option pricing as implemented in fOptions.  The third uses the binomial pricing model from fOptions.  The second and third methods are not implemented here.

The fourth method, "interp", creates a linear interpolation of the cdf of returns, calculates Omega as a vector, and finally interpolates a function for Omega as a function of L.  This method requires library \code{Hmisc}, which can be found on CRAN.
}
\value{
}
\references{
Keating, J. and Shadwick, W.F. The Omega Function. working paper. Finance Development Center, London. 2002.
Kazemi, Schneeweis, and Gupta. Omega as a Performance Measure. 2003.
}
\author{ Peter Carl }
\seealso{
   \code{\link[Hmisc]{Ecdf}}
}
\examples{
    data(edhec)

    for(i in 1:13) { print(Omega(edhec[,i], method="simple")) }

    Omega(edhec[,13],method="interp",output="point")

    Omega(edhec[,13],method="interp",output="full")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ multivariate }
\keyword{ distribution }
\keyword{ models }
