\name{table.Arbitrary}
\alias{table.Arbitrary}
\alias{statsTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ wrapper function for combining arbitrary function list into a table }
\description{
This function creates a table of statistics from vectors of functions and labels passed in.  The resulting table is formatted such that metrics are calculated separately for each column of returns in the data object.

Assumes an input of period returns.  Scale arguements can be used to specify the number of observations during a year (e.g., 12 = monthly returns).
}
\usage{
table.Arbitrary(R, metrics = c("mean", "sd"), metricsNames = c("Average Return", "Standard Deviation"), ...)
statsTable(R, metrics = c("mean", "sd"), metricsNames = c("Average Return", "Standard Deviation"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{ an xts, vector, matrix, data frame, timeSeries or zoo object of asset returns }
  \item{metrics}{ lisdt of functions to apply }
  \item{metricsNames}{ column names for each function }
  \item{\dots}{ any other passthru parameters }
}
\details{
The idea here is to be able to pass in sets of metrics and values, like:

metrics = c(DownsideDeviation(x,MAR=mean(x)), sd(subset(x,x>0)),
sd(subset(x,x<0)), DownsideDeviation(x,MAR=MAR),
DownsideDeviation(x,MAR=Rf=0), DownsideDeviation(x,MAR=0),maxDrawdown(x))

metricsNames = c("Semi Deviation", "Gain Deviation", "Loss Deviation",
paste("Downside Deviation (MAR=",MAR*scale*100,"\%)", sep=""),
paste("Downside Deviation (rf=",rf*scale*100,"\%)", sep=""),
paste("Downside Deviation (0\%)", sep=""), "Maximum Drawdown" )

Here's how it's working right now:
> table.Arbitrary(monthlyReturns.ts,metrics=c("VaR","mean"),
metricsNames=c("modVaR","mean"),p=.95)
\preformatted{
           Actual   S&P500TR
modVaR 0.04186461 0.06261451
mean   0.00945000 0.01013684
}

Passing in two different sets of attributes to the same function doesn't quite work currently.  The issue is apparent in:
> table.Arbitrary(edhec,metrics=c("VaR", "VaR"),
  metricsNames=c("Modified VaR","Traditional VaR"), modified=c(TRUE,FALSE))
\preformatted{
                Convertible.Arbitrage CTA.Global Distressed.Securities
Modified VaR               0.04081599  0.0456767             0.1074683
Traditional VaR            0.04081599  0.0456767             0.1074683
                Emerging.Markets Equity.Market.Neutral Event.Driven
Modified VaR           0.1858624            0.01680917    0.1162714
Traditional VaR        0.1858624            0.01680917    0.1162714
                Fixed.Income.Arbitrage Global.Macro Long.Short.Equity
Modified VaR                 0.2380379   0.03700478        0.04661244
Traditional VaR              0.2380379   0.03700478        0.04661244
                Merger.Arbitrage Relative.Value Short.Selling Funds.of.Funds
Modified VaR          0.07510643     0.04123920     0.1071894     0.04525633
Traditional VaR       0.07510643     0.04123920     0.1071894     0.04525633
}

In the case of this example, you would simply call VaR as the second function, like so:
> table.Arbitrary(edhec,metrics=c("VaR", "VaR"),metricsNames=c("Modified VaR","Traditional VaR"))
\preformatted{
                Convertible.Arbitrage CTA.Global Distressed.Securities
Modified VaR               0.04081599 0.04567670            0.10746831
Traditional VaR            0.02635371 0.04913361            0.03517855
                Emerging.Markets Equity.Market.Neutral Event.Driven
Modified VaR          0.18586240            0.01680917   0.11627142
Traditional VaR       0.07057278            0.01746554   0.03563019
                Fixed.Income.Arbitrage Global.Macro Long.Short.Equity
Modified VaR                0.23803787   0.03700478        0.04661244
Traditional VaR             0.02231236   0.03692096        0.04318713
                Merger.Arbitrage Relative.Value Short.Selling Funds.of.Funds
Modified VaR          0.07510643     0.04123920     0.1071894     0.04525633
Traditional VaR       0.02510709     0.02354012     0.0994635     0.03502065
}
but we don't know of a way to compare the same function side by side with different parameters for each.  Suggestions Welcome.
}
\author{ Peter Carl }
\examples{
data(edhec)
table.Arbitrary(edhec,metrics=c("VaR", "ES"),metricsNames=c("Modified VaR","Modified Expected Shortfall"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ multivariate }
\keyword{ distribution }
\keyword{ models }