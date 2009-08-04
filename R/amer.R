#' fit a (generalized) additive mixed model using lmer
#'
#' @param formula a two-sided formula object describing the fixed-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. The terms can include specifications for non-grouped random effects of types given in \code{basisGenerators}, see examples. The vertical bar character "|" separates an expression for a model matrix and a grouping factor.
#' @param data an optional data frame containing the variables named in formula. By default the variables are taken from the environment from which lmer is called.
#' @param family a GLM family, see \code{\link[stats]{glm}} and \code{\link[stats]{family}}. If \code{family} is missing then a linear mixed model is fit; otherwise a generalized linear mixed model is fit.
#' @param REML logical argument to \code{lmer} only. Should the estimates be chosen to optimize the REML criterion (as opposed to the log-likelihood)?  Defaults to \code{TRUE}.
#' @param control a list of control parameters for \code{\link[lme4]{lmer}}
#' @param start a named list of starting values for the parameters in the model. See \code{\link[lme4]{lmer}}.
#' @param verbose logical scalar.  If \code{TRUE} verbose output is generated during the optimization of the parameter estimates.
#' @param subset see \code{\link[lme4]{lmer}}
#' @param weights see \code{\link[lme4]{lmer}}
#' @param na.action see \code{\link[lme4]{lmer}}
#' @param offset see \code{\link[lme4]{lmer}}
#' @param contrasts see \code{\link[lme4]{lmer}}
#' @param basisGenerators a character vector of names of functions that generate bases for function estimation in  a way lamer can use. See \code{\link{tp}} for an example. 
#'
#' @return An object of class \code{amer-class}.  See there for details.
#' @seealso \code{tests/optionsTests.r} and the vignette for examples 
#' @author Fabian Scheipl
#' @export
amer <-
function(formula, data, family = NULL, REML = TRUE, control = list(), 
		start = NULL, verbose = FALSE, subset, weights, 
		na.action, offset, contrasts = NULL, 
		basisGenerators = c("tp","tpU"),...)# only truncated power basis implemented right now (because tp-bases retain at least some sparsity with the iid penalty)
{	
	expCall <- expand.call()
	
	#set up unfitted model object: 
	setup <- amerSetup(expCall)
	m <- setup$m
	fct <- setup$fct
	fctterm <- setup$fctterm
	rm(setup) #potentially large
	
	#fit the model:
	m <- do.call(
			if(is.null(m$glmFit)){
						lme4:::lmer_finalize
					} else {
						lme4:::glmer_finalize
					}, m)
	
	# add assign-like info to fctterm
	fctterm <- indsF(m, fct, fctterm) 

	ans <- new("amer", m, smooths = fctterm)
	ans@call <- expCall
	# FIXME: see amerSetup: enforce intercept or not?
	###ans@call$formula <- update(formula, .~.+1) #clean up, add intercept to formula
	ans@frame <- merge(ans@frame, data) #add original variables back to frame for prediction etc.
	
	return(ans)
}

