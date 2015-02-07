\name{nwfscDeltaGLM-package}
\alias{nwfscDeltaGLM-package}
\alias{nwfscDeltaGLM}
\docType{package}
\title{
Package for conducting Bayesian estimation of CPUE index standardization (fisheries, ecology, etc).
~~ package title ~~
}
\description{
The nwfscDeltaGLM package is developed by the Northwest Fisheries Science Center, National Oceanic Administration, USA. It's primary purpose is for application to stock assessment, where it can be used to generate time series of CPUE indices from raw trawl data. Optional model components include vessel, strata, and year effects, as well as their interactions. A number of likelihoods are included for the positive component of the model, including both continuous (lognormal, gamma, inverse Gaussian) and discrete (Poisson, Negative Binomial) distributions.
}
\details{
\tabular{ll}{
Package: \tab \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-05\cr
License: \tab What license is it under?\cr
}
}
\author{
Eric J Ward

Maintainer: Eric Ward <eric.ward@noaa.gov>
}
\references{
Thorson, J.T. and E.J. Ward. 2013. Accounting for space-time interactions in index standardization models. Fisheries Research, 147, 426:433. 
}
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%~~ simple examples of the most important functions ~~
#data(Example_Species)
#masterDat = Example_Species
#LoadFn(paste(DataFile,"RandomSpecies.dmp",sep=""))
#strata.limits <- readIn(ncol=5,nlines=6)
#  STRATA  NLat SLat MinDepth MaxDepth
#  A      49.0 36.0  55        183
#  B      49.0 36.0  183       549
#  C      49.0 34.5  549       700
#  D      36.0 34.5  55        549
#  E      34.5 32.0  55        549

# Modify data slightly
#species = "Example_Species"
#names(masterDat)[9] = species

# Preliminary data processing
#processData()

# Define settings
#mcmc.control = list(chains=2, thin=2, burnin=2e3, iterToSave=2e3)
#Parallel = TRUE   # If having trouble, try turning off parallel
#modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="random", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="random", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="fixed", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="fixed", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
#modelStructure2 = list("StrataYear.positiveTows"="correlated", "VesselYear.positiveTows"="correlated", "StrataYear.zeroTows"="correlated", "VesselYear.zeroTows"="correlated", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
#attach(Data)
#mods = list()
#mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1, mcmc.control=mcmc.control,Parallel=Parallel, Species=species)
#mods[[2]] = fitCPUEModel(modelStructure=modelStructure2, mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

# Process MCMC output
# Make sure that Data is attached prior to running
#doMCMCDiags(my.wd,mods)
}
