#' Generate a 2-D truncated power basis
#'
#' Generates a two-dimensional tensor product spline basis of marginal truncated power bases with a single smoothing parameter. 
#' The penalty matrix is the sum of the two kronecker products of the penalties for the marginal bases with the identity matrix.  
#' 
#' @param x covariate for the smooth function
#' @param z covariate for the smooth function
#' @param kx,kz integer: dimensionality of the marginal bases 
#' @param dimU integer: dimension of marginal nullspaces: 1 means unpenalized terms up to linear-linear-interaction, 
#' 			2 means unpenalized terms up to quadratic-quadratic-interactions etc.
#' @param degree integer: degree of (truncated) polynomial
#' @param by factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.
#' @param allPen boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
#'                         if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X
#' @param varying numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying
#' @param diag logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}? 
#' @param knotsx,knotsz double: vector of knot locations (defaults to marginal-quantile-based knot locations)
#' @param centerscale numeric(2): center&scale x by these values if not NULL
#' @param scaledknots boolean:	are knot locations given for the rescaled x-values?
#' @return list with entries:
#'			\code{"X"}: \code{n x (2*degree + degree^2) )} design matrix for unpenalized part (without intercept): marginal polynomials plus their interactions, 
#' 			\code{"Z"}: \code{n x ((kx+1)*(kz+1) - (1 + 2*degree + degree^2))} design matrix for penalized part,
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link{tp}}, \code{\link{tpU}}
tpU2d <- function(x,z, kx=7, kz=7, degree=2, dimU = 1, 
		by=NULL, allPen = FALSE, varying = NULL, diag=FALSE,
		knotsx = quantile(x, probs = (2:(kx - dimU + 1))/(kx - dimU  + 2)),
		knotsz = quantile(z, probs = (2:(kz - dimU + 1))/(kz - dimU  + 2)),
		centerscale=NULL, scaledknots=FALSE)
{
	call <- as.list(expand.call())
	
	#knotsx <- eval(knotsx)
	#knotsz <- eval(knotsz)
	call$knotsx <- knotsx
	call$knotsz <- knotsz
	
	if(is.null(centerscale)){
		x <- scale(x)
		z <- scale(z)
		#make sure prediction data uses the same center/scale as the data used to fit the model:
		call$centerscale <- cbind(x=c(attr(x, "scaled:center"),attr(x, "scaled:scale")),
				z=c(attr(z, "scaled:center"),attr(z, "scaled:scale")))
		x <- as.vector(x)
		z <- as.vector(z)
	} else{
		x <- (x - centerscale[1,1])/centerscale[2,1]
		z <- (z - centerscale[1,2])/centerscale[2,2]
	}	
	if(!scaledknots){
		knotsx <- (knotsx - call$centerscale[1,1])/call$centerscale[2,1]
		knotsz <- (knotsz - call$centerscale[1,2])/call$centerscale[2,2]
		call$scaledknots <- TRUE
	}	
	
	Dx <- cbind(outer(x,0:degree,"^"), outer(x,knotsx,"-")^degree*outer(x,knotsx,">")) 
	
	Dz <- cbind(outer(z,0:degree,"^"), outer(z,knotsz,"-")^degree*outer(z,knotsz,">")) 
	
	##make matrix of tensor product bases (row-wise kronecker-product) 
	#(column-based for efficiency, see mgcv::tensor.prod.model.matrix)
	D <- matrix(0,length(x),0)
	for (j in 1:ncol(Dx)) D<-cbind(D,Dx[,j]*Dz)
	
	
	
	#see also mgcv:::tensor.prod.penalties / DissThomasKneib p. 55
	P <- diag(kz+1)%x%diag(c(rep(0, dimU+1),rep(1,kx-dimU))) + diag(c(rep(0, dimU+1),rep(1,kz-dimU)))%x%diag(kx+1)
	
	res <- reparameterizeDesign(D, P)
	#res$X <- res$X[,-NCOL(res$X), drop=FALSE]#rm intercept
	
	attr(res, "call") <- as.call(call)	
	return(res)	
}

