.onLoad <- function(lib, pkg)
{   
    # Startup Mesage and Desription:
    MSG <- if(getRversion() >= "2.5") packageStartupMessage else message
    dsc <- packageDescription(pkg)
    if(interactive() || getOption("verbose")) { 
        # not in test scripts
        MSG(paste("\nPackage ", pkg, " (",dsc$Version,") loaded.\n",
            # dsc$Title, "\n", 
            "Copyright (c) 2004-2015 Peter Carl and Brian G. Peterson, ", 
            dsc$License, "\n", dsc$URL,
            "\n", sep=""))
    }
}

.onAttach <- function(libname, pkgname) {
  repo <- "https://github.com/braverock/PerformanceAnalytics"
  packageStartupMessage(
    "WARNING: this package was installed from R-Forge, but development has\n",
    "moved to GitHub. Please re-install the package using the GitHub repo at:\n",
    repo, ".")
}

even <- function (x) x%%2==0

odd  <- function (x) x%%2==1

sd.xts <- xts:::sd.xts

#' @importFrom utils packageDescription
#' @importFrom stats sd
#' @import xts
#' @import zoo
NULL

###############################################################################
# R (http://r-project.org/) Econometrics for Performance and Risk Analysis
#
# Copyright (c) 2004-2015 Peter Carl and Brian G. Peterson
#
# This R package is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id$
#
###############################################################################
