#' Reparameterize a spline basis 
#'
#' Reparameterize a spline basis via the spectral decomposition of its penalty matrix into a design for the unpenalized part of the function and 
#' a design for the penalized part with i.i.d. spline coefficients. A centering constraint is applied to the design of the unpenalized part,
#' which only makes sense if the intercept column is actually in the column space of X, i.e. constant functions are in the nullspace of the penalty. 
#' 
#' @param B matrix containing the basis functions
#' @param K penalty matrix (i.e. the penalty term for the fitted spline is the quadratic form of K and the spline coefficients)
#' @param tol eigenvalues smaller than tol are considered zero
#' @return list with entries:
#'			\code{"X"}: \code{n x (p-1)} design matrix for unpenalized part (without intercept) 
#' 			\code{"Z"}: \code{n x (k-p-1)} design matrix for penalized part
#' 			p is the dimension+1 of the penalty nullspace (e.g. 2 for a linear TP-basis or a B-spline with 2nd order difference penalty)  
#' @author Fabian Scheipl
#' @export
#' @seealso source code of \code{tpU2D} for a usage example
reparameterizeDesign <- function(B, K, tol=1e-10)
{
	stopifnot(ncol(B) == ncol(K))
	
	ek <- eigen(K, symmetric = TRUE)
	nullvals <- ek$values < tol
	
	if(any(nullvals)){
		U <- ek$vectors[,nullvals]
		X <- B%*%U
		#absorb centering constraint, see mgcv/smooth.r l.1555 f
		if(ncol(X)>1){
			qrC <- qr(t(matrix(colSums(X),1, ncol(X))))
			X <- t(qr.qy(qrC, t(X)))[,-1, drop=FALSE]
		}else X <- X[, -1, drop=FALSE] 
	} else {
		X <- matrix(0, nr=nrow(B), nc=0)
		warning("penalty matrix seems to have full rank - that's weird!")
	} 
	
	P <- t(1/sqrt(ek$values[!nullvals]) * t(ek$vectors[,!nullvals]))
	Z <- B%*%P
	return(list(X=X, Z=Z))
}	



