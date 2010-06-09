\name{amer}
\alias{amer}
\title{Fit a (generalized) additive mixed model using lmer...}
\usage{amer(formula, data, family=NULL, REML=TRUE, control=list(), start=NULL,
    verbose=FALSE, weights, contrasts=NULL, basisGenerators=c("tp", "tpU",
    "bsp"), ...)}
\description{Fit a (generalized) additive mixed model using lmer}
\value{An object of class \code{\link{amer-class}}.}
\seealso{\code{tests/optionsTests.r} and the vignette for examples}
\author{Fabian Scheipl}
\arguments{\item{formula}{a two-sided formula object describing the fixed-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. The terms can include specifications for non-grouped random effects of types given in \code{basisGenerators}, see examples. The vertical bar character "|" separates an expression for a model matrix and a grouping factor.}
\item{data}{a data frame containing all the variables occuring in the formula.}
\item{family}{a GLM family, see \code{\link[stats]{glm}} and \code{\link[stats]{family}}. If \code{family} is missing then a linear mixed model is fit; otherwise a generalized linear mixed model is fit.}
\item{REML}{logical argument to \code{lmer} only. Should the estimates be chosen to optimize the REML criterion (as opposed to the log-likelihood)?  Defaults to \code{TRUE}.}
\item{control}{a list of control parameters for \code{\link[lme4]{lmer}}}
\item{start}{a named list of starting values for the parameters in the model. See \code{\link[lme4]{lmer}}.}
\item{verbose}{logical scalar.  If \code{TRUE} verbose output is generated during the optimization of the parameter estimates.}
\item{weights}{see \code{\link[lme4]{lmer}}}
\item{contrasts}{see \code{\link[lme4]{lmer}}}
\item{basisGenerators}{a character vector of names of functions that generate bases for function estimation in  a way lamer can use. See \code{\link{tp}} for an example.}
\item{...}{further arguments passed to \code{\link[lme4]{lmer}}.}}
