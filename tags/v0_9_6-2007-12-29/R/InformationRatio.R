`InformationRatio` <-
function (Ra, Rb, scale = 12)
{ # @author Peter Carl

    # DESCRIPTION
    # InformationRatio = ActivePremium/TrackingError

    # Inputs:
    # Outputs:

    # FUNCTION
    assetReturns.vec = checkDataVector(Ra)
    benchmarkReturns.vec = checkDataVector(Rb)

    ActivePremium = ActivePremium(assetReturns.vec,benchmarkReturns.vec, scale = scale)
    TrackingError = TrackingError(assetReturns.vec,benchmarkReturns.vec, scale = scale)

    InformationRatio = ActivePremium/TrackingError

    InformationRatio
}

###############################################################################
# R (http://r-project.org/) Econometrics for Performance and Risk Analysis
#
# Copyright (c) 2004-2007 Peter Carl and Brian G. Peterson
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: InformationRatio.R,v 1.5 2007-04-04 00:23:01 brian Exp $
#
###############################################################################
# $Log: not supported by cvs2svn $
# Revision 1.4  2007/03/11 16:53:19  brian
# - add equations and text to documentation
# - standardize on Ra as the Return of the Asset
# - standardize on Ra as first argument where that wasn't previously true
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