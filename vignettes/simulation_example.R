## ----install-------------------------------------------------------------
# Install package
install.packages("devtools", repos = "http://cran.us.r-project.org")
library("devtools")
install_github("nwfsc-assess/nwfscDeltaGLM")
  
# Load package
library(nwfscDeltaGLM)

## ----load----------------------------------------------------------------
Sim_List = SimSimple(
  MeanEncounter = 0.5,
  Nyears = 10,
  Nstrata = 15,
  Nvessels = 4,
  Obsperyear = 175,
  sigmaV = rep(0, 2),
  sigmaVY = rep(1, 2),
  sigmaS = rep(1, 2),
  sigmaY = rep(1, 2),
  sigmaSY = rep(1, 2),
  sigmaResid = 1
  )

masterDat = Sim_List[["DF"]]

strata.limits = data.frame(
  STRATA = LETTERS[1:15],
  NLat = seq(1.5, 15.5, 1),
  SLat = seq(0.5, 14.5, 1),
  MinDepth = 0,
  MaxDepth = 100,
  stringsAsFactors = FALSE
  )

## ------------------------------------------------------------------------
species = "Simulated_Species"
names(masterDat)[which(names(masterDat)=="SPECIES_WT_KG")] = species
masterDat = cbind(masterDat, AREA_SWEPT_MSQ=masterDat$AREA_SWEPT_HA*1e4)

## ------------------------------------------------------------------------
Covariates = list(positive=TRUE, binomial=TRUE)
nX.binomial = 1
nX.pos = 1
X.bin = matrix( rnorm(nrow(masterDat)*nX.binomial, mean=0, sd=1), ncol=nX.binomial)
X.pos = matrix( rnorm(nrow(masterDat)*nX.pos, mean=0, sd=1), ncol=nX.pos)

## ------------------------------------------------------------------------
DataList = processData()

## ------------------------------------------------------------------------
mcmc.control = list(chains=3, thin=1, burnin=100, iterToSave=200)

modelStructure1 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="random2", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="random2", "Vessel.positiveTows"="zero", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

modelStructure2 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="zero", "Vessel.positiveTows"="random2", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="zero", "Vessel.zeroTows"="random2", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

## ------------------------------------------------------------------------
fitted_models = list()
fitted_models[[1]] = fitDeltaGLM(datalist=DataList, likelihood = "gamma", modelStructure=modelStructure1, mcmc.control=mcmc.control, covariates=Covariates, Species=species)
fitted_models[[2]] = fitDeltaGLM(datalist=DataList, likelihood = "gamma", modelStructure=modelStructure2, mcmc.control=mcmc.control, covariates=Covariates, Species=species)

