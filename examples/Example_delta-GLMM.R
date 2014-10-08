

# Install package
install.packages("devtools")
library("devtools")
install_github("nwfsc-assess/nwfscDeltaGLM")
  
# Load package
library(nwfscDeltaGLM)
# updateDeltaGLMcode()

# File structure
my.wd <- "C:/Users/James.Thorson/Desktop/"
setwd(my.wd)

# Load data and strata
data(Example_Species)
masterDat = Example_Species
strata.limits <- readIn(ncol=5,nlines=6)
  STRATA  NLat SLat MinDepth MaxDepth
  First      49.0 36.0  55        183
  Second      49.0 36.0  183       549
  Mid      49.0 34.5  549       700
  D      36.0 34.5  55        549
  South      34.5 32.0  55        549

# Modify data slightly
species = "Example_Species"
names(masterDat)[9] = species

# Preliminary data processing
processData()

# Define settings
mcmc.control = list(chains=2, thin=1, burnin=1e3, iterToSave=1e3)
Parallel = FALSE   # If having trouble, try turning off parallel
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="randomExpanded", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="randomExpanded", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
mods = list()
mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1, mcmc.control=mcmc.control,Parallel=Parallel, Species=species)

# Process MCMC output
# Make sure that Data is attached prior to running
data(SA3)
doMCMCDiags(my.wd,mods)
