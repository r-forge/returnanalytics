`VaR.CornishFisher` <-
function(R, p=0.99, modified = TRUE)
{   # @author Brian G. Peterson (completed/debugged fn)
    # @author Diethelm Wuertz (prototype function)


    # Description:

    # The limitations of mean Value-at-Risk are well covered in the literature.
    # Laurent Favre and Jose-Antonio Galeano published a paper in the
    # Fall 2002, volume 5 of the Journal of Alternative Investment,
    # "Mean-Modified Value-at-Risk optimization With Hedge Funds",
    # that proposed a modified VaR calculation that takes the higher moments
    # of non-normal distributions (skewness, kurtosis) into account, and
    # collapses to standard (traditional) mean-VaR if the return stream follows a
    # standard distribution.
    # This measure is now widely cited and used in the literature,
    # and is usually referred to as "Modified VaR" or "Modified Cornish-Fisher VaR"

    # Diethelm Wuertz's original function was called monthlyVaR, but did not
    # contain the required modifications to get to a monthly or an annualized number.
    # I have converted it to VaR.CornishFisher, and made the assumption of p=0.99, with an option for p=0.95 and
    # a collapse to normal mean VaR.

    # FUNCTION:

    # compute zc for the probability we want
    if ( p >= 0.51 ) {
        # looks like p was a percent like .99
        p = 1-p
    }
    zc = qnorm(p)

    R = checkData(R, method="matrix")
    columns = ncol(R)
    columnnames=colnames(R)
    # FUNCTION:
    for(column in 1:columns) {
        r = as.vector(na.omit(R[,column]))
        if (!is.numeric(r)) stop("The selected column is not numeric")

        if (modified) {
            s = skewness(r) #skewness of the distribution
            k = kurtosis(r) #(excess) kurtosis
            Zcf = zc + (((zc^2-1)*s)/6) + (((zc^3-3*zc)*k)/24) - (((2*zc^3)-(5*zc)*s^2)/36)
            VaR = mean(r) - (Zcf * sd(r))
            if (eval(VaR<0)){ #eval added to get around Sweave bitching
                warning(c("Cornish-Fisher Expansion produces unreliable result (inverse risk) for column: ",column," : ",VaR))
                # set VaR to 0, since inverse risk is unreasonable
                VaR=0
            }
            if (eval(VaR>1)){ #eval added to get around Sweave bitching
                warning(c("Cornish-Fisher Expansion produces unreliable result (risk over 100%) for column: ",column," : ",VaR))
                # set VaR to 1, since greater than 100% is unreasonable
                VaR=1
            }
        } else {
            VaR = mean(r) - (zc * sd(r))
        }
        VaR=array(VaR)
        if (column==1) {
            #create data.frame
            result=data.frame(VaR=VaR)
        } else {
            VaR=data.frame(VaR=VaR)
            result=cbind(result,VaR)
        }
    } #end columns loop

    if(ncol(result) == 1) {
        # some backflips to name the single column zoo object
        result = as.numeric(result)
    }
    else
        colnames(result) = columnnames

    # Return Value:
    result
}

###############################################################################

`modifiedVaR` <-
function(R, p=0.99)
{   # @author Brian G. Peterson

    # Description:

    # This is a wrapper function for VaR.CornishFisher,
    # because this measure is often referred to as modifiedVaR

    # FUNCTION:
    VaR.CornishFisher(R = R, p = p, modified=TRUE)

}

###############################################################################

`VaR.mean` <-
function(R, p=0.95)
{   # @author Brian G. Peterson

    # Description:

    # This is a wrapper function for modified VaR which assumes a normal
    # distribution by discounting influence from skewness or kurtosis.

    # Wrapper should be used with metrics related to VaR, such as Beyond VaR.

    # FUNCTION:
    VaR.CornishFisher(R = R, p = p, modified=FALSE)

}

###############################################################################

`VaR.traditional` <-
function(R, p=0.95)
{   # @author Brian G. Peterson

    # Description:

    # This is a wrapper function for modified VaR which assumes a normal
    # distribution by discounting influence from skewness or kurtosis.

    # Wrapper should be used with metrics related to VaR, such as Beyond VaR.

    # FUNCTION:
    VaR.CornishFisher(R = R, p = p, modified=FALSE)

}

###############################################################################
# R (http://r-project.org/) Econometrics for Performance and Risk Analysis
#
# Copyright (c) 2004-2007 Peter Carl and Brian G. Peterson
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: VaR.CornishFisher.R,v 1.16 2007-12-29 19:25:09 brian Exp $
#
###############################################################################
# $Log: not supported by cvs2svn $
# Revision 1.14  2007/09/04 02:12:33  brian
# - add eval to if statement for Sweave pickiness
#
# Revision 1.13  2007/07/30 19:06:59  brian
# - fix typo in equation identified by Samantha Kumaran
#
# Revision 1.12  2007/04/04 00:23:01  brian
# - typos and minor comment updates
#
# Revision 1.11  2007/04/02 21:49:22  peter
# - minor modification
#
# Revision 1.10  2007/03/30 14:31:26  peter
# - when a single column is submitted, result is now a "numeric" rather than a
#   "list" object
#
# Revision 1.9  2007/03/22 14:24:15  peter
# - removed column attribute
#
# Revision 1.8  2007/03/22 12:15:25  brian
# - remove VaR.multicolumn, obsolete
#
# Revision 1.7  2007/03/22 11:54:23  brian
# - added handling for multicolumn data
#
# Revision 1.6  2007/03/20 03:26:12  peter
# - removed firstcolumn
#
# Revision 1.5  2007/03/19 21:55:57  peter
# - replaced data checking with checkData function
#
# Revision 1.4  2007/03/11 16:58:07  brian
# - replace as.vector() with checkDataVector()
#
# Revision 1.3  2007/03/04 20:59:27  brian
# - minor changes to pass R CMD check
#
# Revision 1.2  2007/02/07 13:24:49  brian
# - fix pervasive comment typo
#
# Revision 1.1  2007/02/02 19:06:15  brian
# - Initial Revision of packaged files to version control
# Bug 890
#
###############################################################################