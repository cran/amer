\name{amer-class}
\title{Additive Mixed Model class}
\docType{class}
\alias{amer-class}

\description{The \code{"amer"}-class represents additive or generalized
  additive mixed-effects models.  It extends the \code{"mer"}-class.
}

\section{Objects from the Class}{Objects can be created by calls to the
  \code{\link{amer}} function.
}

\section{Slots}{
  The class \code{"amer"} extends \code{\link[lme4]{mer-class}}.
  It has only one additional slot:

  \describe{

    \item{\code{smooths}:}{a ist with one entry for each smooth term in the model. 
    Each entry contains the term in the model formula used to generate the spline basis, with attributes 
    \code{indGrp}, \code{indPen}, \code{indUnpen}, \code{indConst} giving the respective indices of the variance components, 
    random effects, non-constant fixed effects and constant (i.e. intercept, factor level effects)  fixed effects 
    associated with the smooth}
	}
}
\section{Methods}{
	 see \code{\link[lme4]{mer-class}} for a description of methods for \code{mer}-objects
	 
	 \describe{
    \item{predict}{\code{signature(object = "mer")}  See \code{\link[amer]{predict}}}
     } 
}
\seealso{
  \code{\link{getF}} for extracting function values of the smooth terms.
}

\keyword{classes}


