#' returns a call in which all of the arguments which were supplied or have presets are specified by their full names and supplied or default values. 
#' @param definition a function. See \code{\link[base]{match.call}}.
#' @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
#' @param expand.dots logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.
#' @return An object of class call. 
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link[base]{match.call}}
expand.call <-
function(definition=NULL, call=sys.call(sys.parent(1)), expand.dots = TRUE)
{
	call <- .Internal(match.call(definition, call, expand.dots))
	#given args:
	ans <- as.list(call)
	#possible args:
	frmls <- formals(deparse(ans[[1]]))
	#remove formal args with no presets:
	frmls <- frmls[!sapply(frmls, is.symbol)]
	
	add <- which(!(names(frmls) %in% names(ans)))
	return(as.call(c(ans, frmls[add])))
}

