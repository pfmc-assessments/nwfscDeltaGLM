

# Install dependencies
install.packages('rjags')
install.packages('R2jags')
install.packages('pscl')
install.packages('runjags')
install.packages('statmod')
install.packages('superdiag')

# Install package
install.packages("nwfscDeltaGLM", repos="http://R-Forge.R-project.org", type="source")

# Load package
library(nwfscDeltaGLM)
# updateDeltaGLMcode()

# File structure
my.wd <- "C:/Users/James.Thorson/Desktop/"
setwd(my.wd)

# Simulate data
  # MeanEncounter=0.5; Nyears=10; Nstrata=15; Nvessels=4; Obsperyear=175; sigmaV=rep(1,2); sigmaVY=rep(1,2); sigmaS=rep(1,2); sigmaY=rep(1,2); sigmaSY=rep(1,2); sigmaResid=1
Sim_List = SimSimple( MeanEncounter=0.5, Nyears=10, Nstrata=15, Nvessels=4, Obsperyear=175, sigmaV=rep(0,2), sigmaVY=rep(1,2), sigmaS=rep(1,2), sigmaY=rep(1,2), sigmaSY=rep(1,2), sigmaResid=1 )
masterDat = Sim_List[["DF"]]
strata.limits <- readIn(ncol=5,nlines=16)
  STRATA  NLat SLat MinDepth MaxDepth
  A      1.5 0.5  0       100
  B      2.5 1.5  0       100
  C      3.5 2.5  0       100
  D      4.5 3.5  0       100
  E      5.5 4.5  0       100
  F      6.5 5.5  0       100
  G      7.5 6.5  0       100
  H      8.5 7.5  0       100
  I      9.5 8.5  0       100
  J      10.5 9.5  0       100
  K      11.5 10.5  0       100
  L      12.5 11.5  0       100
  M      13.5 12.5  0       100
  N      14.5 13.5  0       100
  O      15.5 14.5  0       100
  
# Generate useless covariates (i.e., not in the true operating model)
Covariates = list(positive=TRUE, binomial=TRUE)
nX.binomial = 1
nX.pos = 1
X.bin = matrix( rnorm(nrow(masterDat)*nX.binomial, mean=0, sd=1), ncol=nX.binomial)
X.pos = matrix( rnorm(nrow(masterDat)*nX.pos, mean=0, sd=1), ncol=nX.pos)

# Modify data slightly
species = "Simulated_Species"
names(masterDat)[5] = species
masterDat = cbind(masterDat, AREA_SWEPT_MSQ=masterDat$AREA_SWEPT_HA*1e4)
processData()

# Define settings
#mcmc.control = list(chains=3, thin=1e1, burnin=1e4, iterToSave=1e4)
mcmc.control = list(chains=2, thin=1e0, burnin=1e3, iterToSave=1e3)
Parallel = FALSE   # If having trouble, try turning off parallel

# Define models
mods = list()
modelStructure1 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="random2", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="random2", "Vessel.positiveTows"="zero", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
  mods[[1]] = fitDeltaGLM(likelihood = "gamma", modelStructure=modelStructure1, mcmc.control=mcmc.control, covariates=Covariates, Parallel=Parallel, Species=species)
modelStructure2 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="zero", "Vessel.positiveTows"="random2", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="zero", "Vessel.zeroTows"="random2", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
  mods[[2]] = fitDeltaGLM(likelihood = "gamma", modelStructure=modelStructure2, mcmc.control=mcmc.control, covariates=Covariates, Parallel=Parallel, Species=species)
save(mods, file=paste(my.wd,"mods.RData",sep=""))

# Process MCMC output
data(SA3)
doMCMCDiags(my.wd, mods, StrataWeights="Equal")
