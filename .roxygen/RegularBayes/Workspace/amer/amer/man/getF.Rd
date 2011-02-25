\name{getF}
\alias{getF}
\title{get the estimated function values from an amer-Fit...}
\usage{getF(object, which, n=100, newdata, interval=c("NONE", "MCMC", "RW"),
    addConst=TRUE, varying=1, level=0.9, sims=1000)}
\description{get the estimated function values from an amer-Fit}
\note{The formula used for the pointwise bias-adjusted CIs is taken from Ruppert and Wand's  'Semiparametric Regression' (2003), p. 140. 
These leave out the uncertainty associated with the variance component estimates. 
MCMC-intervals based on results from \code{\link[lme4]{mcmcsamp}} don't seem to be very reliable yet and should be used with caution, especially for more complex models.}
\value{a list with one \code{data.frame} for each function, giving \code{newdata} or the values of the generated grid plus the fitted values (and confidence/HPD intervals)
if MCMC-intervals were rquested, the list has an attribute "mcmc" containing the result of the call to \code{\link[lme4]{mcmcsamp}}, a \code{\link[lme4]{merMCMC-class}} object.}
\author{Fabian Scheipl}
\seealso{\code{\link{plotF}},  \code{tests/optionsTests.r} and the vignette for examples}
\arguments{\item{object}{a fitted additive (mixed) model of class \code{\link{amer-class}}}
\item{which}{(optional) an integer vector or a character vector of names giving the smooths for which fitted values are desired. Defaults to all.}
\item{n}{if no \code{newdata} is given, fitted values for a regular grid with n values in the range of the respective covariates are returned}
\item{newdata}{An optional data frame in which to look for variables with which to predict}
\item{interval}{what mehod should be used to compute pointwise confidence/HPD intervals: RW= bias-adjusted empirical bayes, MCMC uses \code{\link[lme4]{mcmcsamp}}}
\item{addConst}{boolean should the global intercept and intercepts for the levels of the by-variable be included in the fitted values (and their CIs) can also be a vector of the same length as \code{which}}
\item{varying}{value of the\code{varying}-covariate (see \code{\link{tp}}) to be used if no newdata is supplied. 
Defaults to 1.}
\item{level}{level for the confidence/HPD intervals}
\item{sims}{how many iterates should  be generated for the MCMC-based HPD-intervals}}
