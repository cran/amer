\name{reparameterizeDesign}
\alias{reparameterizeDesign}
\title{Reparameterize a spline basis...}
\usage{reparameterizeDesign(B, K, tol=1e-10)}
\description{Reparameterize a spline basis}
\details{Reparameterize a spline basis via the spectral decomposition of its penalty matrix into a design for the unpenalized part of the function and 
a design for the penalized part with i.i.d. spline coefficients}
\value{list with entries:
\code{"X"}: \code{n x degree} design matrix for unpenalized part (without intercept) (or a list of those for every level of by if allPen=F), 
\code{"Z"}: \code{n x (k-degree)} design matrix for penalized part (or a list of those for every level of by if allPen=F),}
\author{Fabian Scheipl}
\seealso{source code of \code{tp2D} for a usage example}
\arguments{\item{B}{matrix containing the basis functions}
\item{K}{penalty matrix}
\item{tol}{eigenvalues smaller than tol are considered zero}}
