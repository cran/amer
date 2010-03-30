\name{expand.call}
\alias{expand.call}
\title{returns a call in which all of the arguments which were supplied or have presets are specified by their full names and supplied or default values.}
\usage{expand.call(definition, call=sys.call(sys.parent(1)),
    expand.dots=TRUE)}
\description{returns a call in which all of the arguments which were supplied or have presets are specified by their full names and supplied or default values.}
\value{An object of class call.}
\author{Fabian Scheipl}
\seealso{\code{\link[base]{match.call}}}
\arguments{\item{definition}{a function. See \code{\link[base]{match.call}}.}
\item{call}{an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.}
\item{expand.dots}{logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.}}
