## ----install-------------------------------------------------------------
# Install package
install.packages("devtools",repos = "http://cran.us.r-project.org")
library("devtools")
install_github("nwfsc-assess/nwfscDeltaGLM")
  
# Load package
library(nwfscDeltaGLM)

## ----load----------------------------------------------------------------
# Load data and strata
data(Example_Species)
masterDat = Example_Species

strata.limits = data.frame("STRATA" = c("First", "Second", "Mid", "D", "South"),
  NLat = c(49, 49, 49, 36, 34.5), 
  SLat = c(36, 36, 34.5, 34.5, 34.5),
  MinDepth = c(55, 183, 549, 55, 55),
  MaxDepth = c(183, 549, 700, 549, 549), stringsAsFactors = FALSE)

## ------------------------------------------------------------------------
names(masterDat)

species = "Example_Species"
names(masterDat)[which(names(masterDat)=="EXPANDED_WT_KG")] = "Example_Species"

## ------------------------------------------------------------------------
DataList = processData()

## ------------------------------------------------------------------------
mcmc.control = list(chains=1, thin=1, burnin=10, iterToSave=20)

modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="zero", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="zero", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

## ------------------------------------------------------------------------
fitted_models = list()
fitted_models[[1]] = fitDeltaGLM(datalist = DataList, modelStructure=modelStructure1, mcmc.control=mcmc.control, Species=species)

## ----eval=FALSE----------------------------------------------------------
#  # Make sure that Data is attached prior to running
#  data(SA3)
#  doMCMCDiags(datalist=DataList, raw_data = Data, mods=fitted_models, strata.limits=strata.limits)

