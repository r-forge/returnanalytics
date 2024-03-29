`CAPM.beta` <-
function (R, Rb, rf = 0, digits = 4)
{ # @author Peter Carl

    # DESCRIPTION:
    # This is a wrapper for calculating a CAPM beta.

    # Inputs:
    # R: vector of returns for the asset being tested
    # Rb: vector of returns for the benchmark the asset is being gauged against
    # R and Rb are assumed to be matching periods and NO TEST IS MADE TO CHECK
    # rf: risk free rate in the same periodicity as the returns.  May be a vector
    #     of the same length as x and y.

    # Output:
    #

    # FUNCTION:

    assetReturns.vec = checkDataVector(R)
    benchmarkReturns.vec = checkDataVector(Rb)
    rf.vec = checkDataVector(rf)

    if (length(assetReturns.vec) != length(benchmarkReturns.vec))
        stop("Returns to be assessed have unequal time periods. Are there NA's in the data?")

    # Make these excess returns
    assetExcessRet.vec = assetReturns.vec - rf.vec
    indexExcessRet.vec = benchmarkReturns.vec - rf.vec

    # regress
    model.lm = lm(assetExcessRet.vec ~ indexExcessRet.vec)

#    alpha = coef(model.lm)[[1]]
    beta = coef(model.lm)[[2]]
    beta
}

###############################################################################
# R (http://r-project.org/) Econometrics for Performance and Risk Analysis
#
# Copyright (c) 2004-2007 Peter Carl and Brian G. Peterson
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: CAPM.beta.R,v 1.5 2007-02-28 03:26:00 peter Exp $
#
###############################################################################
# $Log: not supported by cvs2svn $
# Revision 1.4  2007/02/08 21:43:39  brian
# - standardize parameters to R and Rb for consistency with other functions
#
# Revision 1.3  2007/02/07 13:24:49  brian
# - fix pervasive comment typo
#
# Revision 1.2  2007/02/07 13:20:52  brian
# - change Ri to Rb for benchmark asset to standardize parameters
# - change indexReturns.vec to benchmarkReturns.vec for consistency
#
# Revision 1.1  2007/02/02 19:06:15  brian
# - Initial Revision of packaged files to version control
# Bug 890
#
###############################################################################