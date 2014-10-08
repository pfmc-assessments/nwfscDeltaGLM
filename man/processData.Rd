\name{processData}
\alias{processData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{processData is the function that converts the raw trawl data file into a format that is used by the JAGS model}
\description{processData attaches the appropriate variables to the global workspace}
\usage{
processData(Truncate=0,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Truncate}{An optional truncation value, defaults to 0. Any values below the truncation value will be set = 0}
  \item{\dots}{Other arguments can be passed.}  
}
\value{}
\note{
Two output files are also created with each function call. The file "Tows_outside_strata.csv" is written to the working directory, containing all tows that fall outside the user-definied strata. The file "Tows_with_missing_data.csv" is created, containing all tows that have missing necessary information.
}
\seealso{}
