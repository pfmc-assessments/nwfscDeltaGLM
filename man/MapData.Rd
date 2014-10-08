\name{MapData}
\alias{MapData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MapData is function that makes a map by latitude / depth}
\description{This function produces a plot called TowMap.jpg (originally designed by John Wallace), showing tows in strata by depth and latitude.
}
\usage{
MapData(Data, SA3, strata.limits, FileName, Folder,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{The data frame storing raw tow data}
  \item{SA3}{A flat data file, with area ID, min / max depth, min / max latitude, and the area of each strata in hectares.}
  \item{strata.limits}{A file (e.g. defaultLimits.csv), indicating an alphabetic character for each strata (A, B, C, ...), as well as the depth and latitude ranges for each stratum.}
  \item{FileName}{The name of the output file -- this will be species and model specific}
  \item{Folder}{Where the output is stored}  
  \item{\dots}{Other arguments can be passed.}  
}
\value{}
\seealso{}
