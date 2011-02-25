\name{predict}
\docType{methods}
\alias{predict,amer-method}
\title{Get predicted/fitted values for new data}
\usage{
\S4method{predict}{amer}(object, newdata, type=c("response", "linpred", "terms"), ...)
}
\arguments{
  \item{object}{a fitted model object inheriting from
    class \code{\linkS4class{amer}}.
  }
  \item{newdata}{a \code{data.frame} containing observations with all
   variables used in the original fit for which predicted/fitted values should be computed.}
  \item{type}{see Value}
  \item{\dots}{Additional, optional arguments for some methods.  At
    present none are used.}
}
\description{
 Get predicted/fitted values for the new data supplied. 
}
\value{
  the \code{data.frame} supplied in \code{newdata} with additional columns containing:
  
  the components of the linear predictor \eqn{X\beta}, named as in \code{fixef(object)} (prefixed by \code{lp.}),
  
  the components of the random effects \eqn{Z_i b_i}, named after the grouping factors (prefixed by \code{blup.}),
  
  the smooth terms, named as in \code{names(objects@smooths)}.
  
  If \code{type="response"} the \code{"fit"} column contains the fitted values on the scale of the response,
  if \code{type="linpred"} the \code{"fit"} column contains the linear predictor (i.e. fixed + random + smooth effects).
}  
\note{ 
  If \code{newdata} contains additional levels of grouping factors not present in the original data, 
  these are assigned random effects of zero.
  
  The function doesn't check whether the new data supplied for the smooth terms is within the range of the original data. 
  Extrapolating spline fits can be dangerous.
}
\author{Fabian Scheipl}

