\name{nwfscDeltaGLM-package}
\alias{nwfscDeltaGLM-package}
\alias{nwfscDeltaGLM}
\docType{package}
\title{
Package for conducting Bayesian estimation of CPUE index standardization (fisheries, ecology).
}
\description{
The nwfscDeltaGLM package is developed by the Northwest Fisheries Science Center, National Oceanic Administration, USA. It's primary purpose is for application to stock assessment, where it can be used to generate time series of CPUE indices from raw trawl data. Optional model components include vessel, strata, and year effects, as well as their interactions. A number of likelihoods are included for the positive component of the model, including both continuous (lognormal, gamma, inverse Gaussian) and discrete (Poisson, Negative Binomial) distributions.
}
\author{
Maintainer: Eric Ward <eric.ward@noaa.gov>
}
\references{
Thorson, J.T. and E.J. Ward. 2013. Accounting for space-time interactions in index standardization models. Fisheries Research, 147, 426-433.

Thorson, J.T., A.O. Shelton, E.J. Ward, and H.J. Skaug. 2015. Geostatistical delta-generalized linear mixed models improve precision for estimated abundance indices for West Coast groundfishes. ICES J Mar Sci, 72, 1297-1310.
}
\keyword{ package }
\examples{

data(Example_Species) # raw haul data
data(SA3) # raw strata areas
data(defaultLimits) # default limits

masterDat = Example_Species

# Modify data slightly
species = "Example_Species"
names(masterDat)[which(names(masterDat)=="EXPANDED_WT_KG")] = species

# Preliminary data processing
DataList = processData(masterDat = masterDat, strata.limits=defaultLimits)

# Define settings
mcmc.control = list(chains=3, thin=1, burnin=500, iterToSave=1000)

modelStructure1 = list("StrataYear.positiveTows"="fixed",
"VesselYear.positiveTows"="zero",
"StrataYear.zeroTows"="fixed",
"VesselYear.zeroTows"="zero",
"Vessel.positiveTows"="zero",
"Vessel.zeroTows"="zero",
"Catchability.positiveTows"="one",
"Catchability.zeroTows"="zero",
"year.deviations"="fixed",
"strata.deviations"="fixed")

# Define models
fitted_models = list()
fitted_models[[1]] = fitDeltaGLM(datalist=DataList, modelStructure=modelStructure1,
mcmc.control=mcmc.control,
Species=species,
Parallel=FALSE)

#mods[[2]] = fitCPUEModel(modelStructure=modelStructure2,
#mcmc.control=mcmc.control,
#Parallel=Parallel,
#Species=species)

# Process MCMC output
# Make sure that Data is attached prior to running
#doMCMCDiags(my.wd,mods)
}
