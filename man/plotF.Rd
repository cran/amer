\name{plotF}
\alias{plotF}
\title{plot estimated function values + pointwise confidence intervals for an amer-Fit...}
\usage{plotF(object, which, n=100, interval="RW", addConst=TRUE,
    trans=I, level=0.9, sims=1000, auto.layout=TRUE, rug=TRUE,
    legendPos="topright", ...)}
\description{plot estimated function values + pointwise confidence intervals for an amer-Fit}
\value{invisibly returns the result of the \code{\link{getF}}-call used to do the plots}
\author{Fabian Scheipl}
\seealso{\code{\link{getF}}; \code{\link{plotF}},  \code{tests/optionsTests.r} and the vignette for examples}
\arguments{\item{object}{a fitted additive (mixed) model of class \code{\link{amer-class}}}
\item{which}{(optional) an integer vector or a character vector of names giving the smooths for which plots are desired. Defaults to all.}
\item{n}{fitted values for a regular grid with n values in the range of the respective covariates are calculated}
\item{interval}{what mehod should be used to compute pointwise confidence/HPD intervals. See \code{\link{getF}}.}
\item{addConst}{boolean should the global intercept and intercepts for the levels of the by-variable be included in the fitted values (and their CIs)}
\item{trans}{a function that should be applied to the fitted values and ci's before plotting (e.g. the inverse link function to get plots on the scale of the reponse)}
\item{level}{level for the confidence/HPD intervals}
\item{sims}{how many iterates should  be generated for the MCMC-based HPD-intervals}
\item{auto.layout}{automagically set plot layout via \code{par()$mfrow}}
\item{rug}{add \code{\link{rug}}-plots of the observed covariate locations}
\item{legendPos}{a (vector of) keyword(s) where to put labels of by-variables (see \code{\link[graphics]{legend}}). "none" if you don't want a legend.}
\item{...}{arguments passed on to the low-level plot functions (\code{plot}, \code{matlines}), \code{legend}, and \code{title}}}
