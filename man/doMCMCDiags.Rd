\name{doMCMCDiags}
\alias{doMCMCDiags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{doMCMCDiags is the master diagnostic function}
\description{This function is a high level function to be called by the user, and implements all of the plotting, diagnostic, and output summaries available. It is applicable to one or more fitted model objects, stored in a list object.
}
\usage{
doMCMCDiags(datalist, strata.limits=strata.limits, directory, mods, StrataWeights="StrataAreas", McmcDiagnostics=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directory}{The directory to save output in}
  \item{mods}{A list of one or more fitted delta-GLMM model objects}
  \item{McmcDiagnostics}{Boolean indicating whether or not to make the MCMCdiagnostics plots, defaults to FALSE.}
  \item{\dots}{Other arguments can be passed.}  
}
\value{}
\seealso{\code{\link{MapData}}, \code{\link{ConvergencePlot}}, \code{\link{McmcDiagnosticsPlot}}, \code{\link{PlotOffset}}, \code{\link{PosteriorPredictive}}, \code{\link{ComputeIndices}}, \code{\link{ComputeMleIndices}}}
