\name{ComputeMleIndices}
\alias{ComputeMleIndices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes maximum likelihood indices from fitted delta-GLMM model}
\description{Takes data from a population census, creates a table of animals only alive at each time period, then fits a Bayesan state space model to the data. By default, a plot of population counts + estimated states is returned, but plots of posteriors of interest (including trend and carrying capacity) can also be generated.}
\usage{
ComputeMleIndices(Data, Model, FileName, Folder=NA, Weights="StrataAreas", StrataTable, Run=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{}
  \item{Model}{}
  \item{FileName}{}
  \item{Folder}{}
  \item{Weights}{}  
  \item{StrataTable}{}  
  \item{Run}{}    
  \item{\dots}{Other arguments can be passed.}  
}
\value{}
\seealso{}
