

# Install dependencies
install.packages('rjags')
install.packages('R2jags')
install.packages('pscl')
install.packages('runjags')
install.packages('statmod')
install.packages('superdiag')

# Install package
install.packages("nwfscDeltaGLM", repos="http://R-Forge.R-project.org")

# Load package
library(nwfscDeltaGLM)
# updateDeltaGLMcode()

# File structure
my.wd <- "C:/Users/James.Thorson/Desktop/"
setwd(my.wd)

# Settings
Nstrata = 15
Nyears = 10
Nrep = 100

# Save object
IndexResults = array(NA, dim=c(Nrep,Nyears,3), dimnames=list( paste("Rep",1:Nrep), paste("Year",1:Nyears), c("True","IndexMedian","SdLog") ))
# load(file="IndexResults.RData")

# Define strata
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

# Loop
for(RepI in 1:Nrep){

  # Simulate data
  Sim_List = SimSimple( MeanEncounter=0.5, Nyears=Nyears, Nstrata=Nstrata, Nvessels=4, Obsperyear=175, sigmaV=rep(0,2), sigmaVY=rep(1,2), sigmaS=rep(1,2), sigmaY=rep(1,2), sigmaSY=rep(1,2), sigmaResid=1, nX.pos=nX.pos, nX.binomial=nX.binomial )
  masterDat = Sim_List[["DF"]]
    
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
  modelStructure1 = list("StrataYear.positiveTows"="random2", "VesselYear.positiveTows"="random2", "StrataYear.zeroTows"="random2", "VesselYear.zeroTows"="random2", "Vessel.positiveTows"="zero", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
    mods[[1]] = fitDeltaGLM(likelihood = "gamma", modelStructure=modelStructure1, mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
  
  # Process MCMC output
  data(SA3)
  doMCMCDiags(my.wd, mods, StrataWeights="Equal")

  # Estimated
  Index = read.csv(paste("Simulated_Species_FinalDiagnostics/Model=1/ResultsByYear.csv",sep=""), header=TRUE)
  IndexResults[RepI,,c("IndexMedian","SdLog")] = as.matrix(Index[,c("IndexMedian","SdLog")])

  # True
  Pres_True = plogis( outer(rep(1,Nstrata),Sim_List$betaY[1,]) + outer(Sim_List$betaS[1,],rep(1,Nyears)) + Sim_List$betaSY[1,,] ) 
  Pos_True = exp( outer(rep(1,Nstrata),Sim_List$betaY[2,]) + outer(Sim_List$betaS[2,],rep(1,Nyears)) + Sim_List$betaSY[2,,] ) 
  Index_True = colSums( Pres_True * Pos_True )
  IndexResults[RepI,,"True"] = Index_True
}

# Credible interval coverage
png( file="CIC.png", width=4, height=4, res=200, units="in")
  hist( (log(IndexResults[,,'IndexMedian'])-log(IndexResults[,,'True']))/IndexResults[,,'SdLog'], freq=FALSE)
  lines( x=seq(-5,5,length=1000), dnorm(seq(-5,5,length=1000)))
dev.off()
