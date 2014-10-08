\name{McmcDiagnosticsPlot}
\alias{McmcDiagnosticsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{McmcDiagnosticsPlot does Gelman-Rubin and optionally Geweke diagnostic plots}
\description{This function does MCMC diagnostics, producing several plots related to Gelman diagnostics, and optionally a plot for the Geweke statistics (though this doesn't converge sometimes).
}
\usage{
McmcDiagnosticsPlot(McmcList, parToMonitor, FileName, Folder=NA, Geweke=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{McmcList}{An MCMC list object of model parameters}
  \item{parToMonitor}{Which parameter(s) of the fitted model object were monitored / estimated}
  \item{FileName}{Species and model specific, where the outupt will be generated to}
  \item{Folder}{Folder for output}  
  \item{Geweke}{Whether to also calculate Geweke diagnostics, defaults to FALSE}
  \item{\dots}{Other arguments can be passed.}  
}
\value{}
\seealso{}
