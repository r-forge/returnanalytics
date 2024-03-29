`sd.multiperiod` <-
function (x, na.rm=FALSE, periods = 12, ...)
{
    x = checkData (x,na.rm=na.rm, ...=...)
    #scale standard deviation by multiplying by the square root of the number of periods to scale by
	return(sqrt(periods)*sd(x, na.rm=na.rm))
}

`sd.annualized` <-
function (x, na.rm=FALSE, periods = 12, ...)
{
	sd.multiperiod(x, na.rm=na.rm, periods = periods, ...=...)
}

`StdDev.annualized` <-
function (Ra, na.rm=FALSE, scale = 12, ...)
{   # wrapper function for backwards compatibility
	sd.multiperiod(Ra, na.rm=na.rm, periods = scale, ...=...)
}

###############################################################################
# sd function wrappers for backwards compatibility
`StdDev` <-
function(Ra)
{ # wrapper for backwards compatibility
    return(sd(Ra))
}

`std` <-
function(Ra) {
    # NOTE: std function is listed in the doc for fBasics, but not implemented
    return(sd(Ra))
}

###############################################################################
# R (http://r-project.org/) Econometrics for Performance and Risk Analysis
#
# Copyright (c) 2004-2007 Peter Carl and Brian G. Peterson
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: StdDev.annualized.R,v 1.8 2007-06-07 23:02:20 brian Exp $
#
###############################################################################
# $Log: not supported by cvs2svn $
# Revision 1.7  2007/06/04 14:32:33  peter
# - fixed x and Ra replacement error
#
# Revision 1.6  2007/05/15 20:02:01  brian
# - fix syntax error (extra paren)
#
# Revision 1.5  2007/05/15 11:57:52  brian
# - standardize usage to match common R usage
# - define sd.annualized and sd.multiperiod as generic extensions of R core sd fn
# - move StdDev and std wrappers to this file
#
# Revision 1.4  2007/03/12 15:45:50  brian
# - add equations to documentation
# - standardize on Ra for Returns of asset
#
# Revision 1.3  2007/02/15 01:14:43  brian
# - standardize parameter variaable names
#
# Revision 1.2  2007/02/07 13:24:49  brian
# - fix pervasive comment typo
#
# Revision 1.1  2007/02/02 19:06:15  brian
# - Initial Revision of packaged files to version control
# Bug 890
#
###############################################################################
