fitCPUEModel = function(modelStructure = list("StrataYear.positiveTows" = "random","VesselYear.positiveTows" = "random","StrataYear.zeroTows" ="random","VesselYear.zeroTows" = "random", "Vessel.positiveTows"="zero", "Vessel.zeroTows"="zero", "Catchability.positiveTows" = "one", "Catchability.zeroTows" = "zero", "year.deviations" = "uncorrelated","strata.deviations" = "uncorrelated"),covariates=list(positive=FALSE,binomial=FALSE),likelihood = "gamma", pres_link="logit", model.name = "deltaGLM.txt", fit.model=TRUE, write.model=TRUE, mcmc.control = list(chains = 5, thin = 1, burn = 5000, iterToSave = 2000),Parallel=TRUE, Species = "NULL",logitBounds = c(-20,20),logBounds = c(-20,20), prior.scale = rep(25,6), dgammaNum=0.001) {

  if(modelStructure$Catchability.positiveTows%in%c("linear","quadratic") | modelStructure$Catchability.zeroTows%in%c("one","linear","quadratic")){
    print("Warning: index will not have comparable scale to a design-based (raw) index unless catchability.positiveTows equals 'one'  catchability.zeroTows equals 'zero'") 
  }
  
  if(.Platform$OS.type != "windows" & Parallel == TRUE) {
  	print("Warning: system OS is non-windows, and running in parallel is currently not supported. Code will run in non-parallel")
  	Parallel = FALSE
  }
  if(length(prior.scale) != 6 | length(which(is.na(prior.scale)==T)) > 0 | length(which(prior.scale <= 0)) > 0) {
  	print("Error: prior.scale needs to be specified as a 6-element vector, of positive values")
  	stop()
  }
  if(Species == "NULL") {
  	print("Error: you need to input a species name in the function, e.g. Species = \"Canary\"")
  	stop()
  }
  if(modelStructure$VesselYear.zeroTows =="correlated" & modelStructure$VesselYear.positiveTows != "correlated") {
  	print("Error: you specified only Vessel Year deviations for binomial model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$VesselYear.zeroTows !="correlated" & modelStructure$VesselYear.positiveTows == "correlated") {
  	print("Error: you specified only Vessel Year deviations for positive model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$Vessel.zeroTows =="correlated" & modelStructure$Vessel.positiveTows != "correlated") {
  	print("Error: you specified only Vessel deviations for binomial model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$Vessel.zeroTows !="correlated" & modelStructure$Vessel.positiveTows == "correlated") {
  	print("Error: you specified only Vessel deviations for positive model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }  
  if(modelStructure$StrataYear.zeroTows =="correlated" & modelStructure$StrataYear.positiveTows != "correlated") {
  	print("Error: you specified only Strata Year deviations for binomial model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$StrataYear.zeroTows !="correlated" & modelStructure$StrataYear.positiveTows == "correlated") {
  	print("Error: you specified only Strata Year deviations for positive model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  
  if(nVY == nY) {
  	# EW: this was a catch I put in for the darkblotched data on 5.13.12. 
  	# When there's only 1 vessel, and VY interaction needs to be set to 0
  	# because it contains redundant info as year
  	modelStructure$VesselYear.positiveTows == "zero"
  	modelStructure$VesselYear.zeroTows == "zero"
  }
  
  # check to make sure that the covariates are in matrix form - this is just a BUGS/JAGS thing
  nX.binomial = 1
  if(covariates$binomial) {
  	# check for matrix
  	if(is.matrix(X.bin)==FALSE) {
  		print("Error: Please specify the covariates for the binomial model as a matrix, or turn covariates off.")
  	  stop()
  	}
  	nX.binomial = dim(X.bin)[2]
  }
  nX.pos = 1
  if(covariates$positive) {
  	# check for matrix
  	if(is.matrix(X.pos)==FALSE) {
  		print("Error: Please specify the covariates for positive tows as a matrix, or turn covariates off.")
  	  stop()
  	}
  	nX.pos = dim(X.pos)[2]	
  }
  
  ########################################################################################################################################
  # These 6 strings are all new, relevant for implementing the variance manipulated/expanded method described in Gelman et al. 2007, 
  # Gelman et al. 2006. The prior implicit on random effects is the non-central folded t distribution from Gelman et al. 2006
  ########################################################################################################################################
  SYexpanded = paste("   tau.xi[1] <- pow(",prior.scale[1],", -2);\n","   xi[1] ~ dnorm (0, tau.xi[1]);\n",
  "tau.eta[1] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaSY[1] <- abs(xi[1])/sqrt(tau.eta[1]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nSY){\n",
  "   SYeta[j] ~ dnorm(0, tau.eta[1]); # hierarchical model for theta\n",
  "   SYdev[j] <- min(max(xi[1]*SYeta[j],",logBounds[1],"),",logBounds[2],");\n",
  "}\n",sep="")
  pSYexpanded = paste("   tau.xi[2] <- pow(",prior.scale[2],", -2);\n","   xi[2] ~ dnorm (0, tau.xi[2]);\n",
  "tau.eta[2] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaSY[2] <- abs(xi[2])/sqrt(tau.eta[2]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nSY){\n",
  "   pSYeta[j] ~ dnorm(0, tau.eta[2]); # hierarchical model for theta\n",
  "   pSYdev[j] <- min(max(xi[2]*pSYeta[j],",logitBounds[1],"),",logitBounds[2],");\n",
  "}\n",sep="")
  VYexpanded = paste("   tau.xi[3] <- pow(",prior.scale[3],", -2);\n","   xi[3] ~ dnorm (0, tau.xi[3]);\n",
  "tau.eta[3] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaVY[1] <- abs(xi[3])/sqrt(tau.eta[3]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nVY){\n",
  "   VYeta[j] ~ dnorm(0, tau.eta[3]); # hierarchical model for theta\n",
  "   VYdev[j] <- min(max(xi[3]*VYeta[j],",logBounds[1],"),",logBounds[2],");\n",
  "}\n",sep="")
  pVYexpanded = paste("   tau.xi[4] <- pow(",prior.scale[4],", -2);\n","   xi[4] ~ dnorm (0, tau.xi[4]);\n",
  "tau.eta[4] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaVY[2] <- abs(xi[4])/sqrt(tau.eta[4]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nVY){\n",
  "   pVYeta[j] ~ dnorm(0, tau.eta[4]); # hierarchical model for theta\n",
  "   pVYdev[j] <- min(max(xi[4]*pVYeta[j],",logitBounds[1],"),",logitBounds[2],");\n",
  "}\n",sep="")
  Vexpanded = paste("   tau.xi[5] <- pow(",prior.scale[5],", -2);\n","   xi[5] ~ dnorm (0, tau.xi[5]);\n",
  "tau.eta[5] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaV[1] <- abs(xi[5])/sqrt(tau.eta[5]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nV){\n",
  "   Veta[j] ~ dnorm(0, tau.eta[5]); # hierarchical model for theta\n",
  "   Vdev[j] <- min(max(xi[5]*Veta[j],",logBounds[1],"),",logBounds[2],");\n",
  "}\n",sep="")
  pVexpanded = paste("   tau.xi[6] <- pow(",prior.scale[6],", -2);\n","   xi[6] ~ dnorm (0, tau.xi[6]);\n",
  "tau.eta[6] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaV[2] <- abs(xi[6])/sqrt(tau.eta[6]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nV){\n",
  "   pVeta[j] ~ dnorm(0, tau.eta[6]); # hierarchical model for theta\n",
  "   pVdev[j] <- min(max(xi[6]*pVeta[j],",logitBounds[1],"),",logitBounds[2],");\n",
  "}\n",sep="")

  ########################################################################################################################################
  # This section is related to strata-year and vessel-year effects
  # Strata-year interactions can be (1) fixed, (2) random, (3) randomExpanded, or (4) not estimated (set to 0)
  ########################################################################################################################################  
  SYpos.string = ""
  if(modelStructure$StrataYear.positiveTows == "fixed") SYpos.string = paste("   for(i in 1:nSY) {\n      SYdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n   }\n   strataYearTau[1,1] <- 0;\n","   strataYearTau[1,2] <- 0;\n   sigmaSY[1]<-0;\n   tauSY[1]<-0;\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "random") SYpos.string = paste("   for(i in 1:nSY) {\n      SYdev[i] ~ dnorm(0,tauSY[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   strataYearTau[1,1] <- 0;\n","   strataYearTau[1,2] <- 0;\n   sigmaSY[1]~dunif(0,100);\n   tauSY[1]<-pow(sigmaSY[1],-2);\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "randomExpanded") SYpos.string = paste(SYexpanded, "   strataYearTau[1,1] <- 0;\n","   strataYearTau[1,2] <- 0;\n   tauSY[1]<-pow(sigmaSY[1],-2);\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "zero") SYpos.string = "   for(i in 1:nSY) {\n      SYdev[i] <- 0;\n   }\n   strataYearTau[1,1] <- 0;\n   strataYearTau[1,2] <- 0;\n   sigmaSY[1]<-0;\n   tauSY[1]<-0;\n"
  
  SYzero.string = ""
  if(modelStructure$StrataYear.zeroTows == "fixed") SYzero.string = paste("   for(i in 1:nSY) {\n      pSYdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n   }\n   strataYearTau[2,1] <- 0;\n","   strataYearTau[2,2] <- 0;\n","   sigmaSY[2] <- 0;\n   tauSY[2]<-0;\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "random") SYzero.string = paste("   for(i in 1:nSY) {\n      pSYdev[i] ~ dnorm(0,tauSY[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   strataYearTau[2,1] <- 0;\n","   strataYearTau[2,2] <- 0;\n   sigmaSY[2]~dunif(0,100);\n   tauSY[2]<-pow(sigmaSY[2],-2);\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "randomExpanded") SYzero.string = paste(pSYexpanded, "   strataYearTau[2,1] <- 0;\n","   strataYearTau[2,2] <- 0;\n   tauSY[2]<-pow(sigmaSY[2],-2);\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "zero") SYzero.string = "   for(i in 1:nSY) {\n      pSYdev[i] <- 0;\n   }\n   strataYearTau[2,1] <- 0;\n   strataYearTau[2,2] <- 0;\n   sigmaSY[2] <- 0;\n   tauSY[2]<-0;\n"
  
  # combine the strata year interactions into a string
  stratayear.string = paste(SYpos.string,SYzero.string)
  if(modelStructure$StrataYear.zeroTows == "correlated" & modelStructure$StrataYear.positiveTows == "correlated") {
  	# strata year deviations are MVN RE
  	stratayear.string = paste("sigmaSY[1] <- 0;\n   sigmaSY[2] <- 0;\n      strataYearTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nSY) {\n   sydevs[i,1:2] ~ dmnorm(zs[1:2],strataYearTau[1:2,1:2]);\n      SYdev[i] <- min(max(sydevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pSYdev[i] <- min(max(sydevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")	
  }
  
  # Vessel-year interactions cab be (1) fixed, (2) random, (3) randomExpanded, or (4) not estimated (set to 0) 
  VYpos.string = ""
  if(modelStructure$VesselYear.positiveTows == "fixed") VYpos.string = paste("   VYdev[1] <- 0;\n   for(i in 2:nVY) {\n      VYdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n   }\n   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]<-0;\n   tauVY[1]<-0;\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "random") VYpos.string = paste("   for(i in 1:nVY) {\n      VYdev[i] ~ dnorm(0,tauVY[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]~dunif(0,100);\n   tauVY[1]<-pow(sigmaVY[1],-2);\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "random2") VYpos.string = paste("   for(i in 1:nVY) {\n      VYdev[i] ~ dnorm(0,tauVY[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]<-pow(tauVY[1],-1/2);\n   tauVY[1]~dgamma(",dgammaNum,",",dgammaNum,");\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "randomExpanded") VYpos.string = paste(VYexpanded,"   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   tauVY[1]<-pow(sigmaVY[1],-2);\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "zero") VYpos.string = "   for(i in 1:nVY) {\n      VYdev[i] <- 0;\n   }\n   vesselYearTau[1,1] <- 0;\n   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]<-0;\n   tauVY[1]<-0;\n"
  
  VYzero.string = ""
  if(modelStructure$VesselYear.zeroTows == "fixed") VYzero.string = paste("   pVYdev[1] <- 0;\n   for(i in 2:nVY) {\n      pVYdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   sigmaVY[2]<-0;\n   tauVY[2]<-0;\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "random") VYzero.string = paste("   for(i in 1:nVY) {\n      pVYdev[i] ~ dnorm(0,tauVY[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   sigmaVY[2] ~ dunif(0,100);\n   tauVY[2]<-pow(sigmaVY[2],-2);\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "random2") VYzero.string = paste("   for(i in 1:nVY) {\n      pVYdev[i] ~ dnorm(0,tauVY[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   sigmaVY[2]<-pow(tauVY[2],-1/2);\n   tauVY[2]~dgamma(",dgammaNum,",",dgammaNum,");\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "randomExpanded") VYzero.string = paste(pVYexpanded,"   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   tauVY[2]<-pow(sigmaVY[2],-2);\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "zero") VYzero.string = "   for(i in 1:nVY) {\n      pVYdev[i] <- 0;\n   }\n   vesselYearTau[2,1] <- 0;\n   vesselYearTau[2,2] <- 0;\n   sigmaVY[2]<-0;\n   tauVY[2]<-0;\n"
  
  # combine the strata year interactions into a string
  vesselyear.string = paste(VYpos.string,VYzero.string)

  if(modelStructure$VesselYear.zeroTows == "correlated" & modelStructure$VesselYear.positiveTows == "correlated") {
  	# vessel year deviations are MVN RE
  	vesselyear.string = paste("sigmaVY[1] <- 0;\n   sigmaVY[2] <- 0;\n   vesselYearTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nVY) {\n   vydevs[i,1:2] ~ dmnorm(zs[1:2],vesselYearTau[1:2,1:2]);\n      VYdev[i] <- min(max(vydevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pVYdev[i] <- min(max(vydevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")	
  }
  
  # Vessel interactions cab be (1) fixed, (2) random, (3) randomExpanded, or (4) not estimated (set to 0) 
  Vpos.string = ""
  if(modelStructure$Vessel.positiveTows == "fixed") Vpos.string = paste("   Vdev[1] <- 0;\n   for(i in 2:nV) {\n      Vdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n   }\n   vesselTau[1,1] <- 0;\n","   vesselTau[1,2] <- 0;\n   sigmaV[1]<-0;\n   tauV[1]<-0;\n",sep="")
  if(modelStructure$Vessel.positiveTows == "random") Vpos.string = paste("   for(i in 1:nV) {\n      Vdev[i] ~ dnorm(0,tauV[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   vesselTau[1,1] <- 0;\n","   vesselTau[1,2] <- 0;\n   sigmaV[1]~dunif(0,100);\n   tauV[1]<-pow(sigmaV[1],-2);\n",sep="")
  if(modelStructure$Vessel.positiveTows == "random2") Vpos.string = paste("   for(i in 1:nV) {\n      Vdev[i] ~ dnorm(0,tauV[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   vesselTau[1,1] <- 0;\n","   vesselTau[1,2] <- 0;\n   sigmaV[1]<-pow(tauV[1],-1/2);\n   tauV[1]~dgamma(",dgammaNum,",",dgammaNum,");\n",sep="")
  if(modelStructure$Vessel.positiveTows == "randomExpanded") Vpos.string = paste(Vexpanded,"   vesselTau[1,1] <- 0;\n","   vesselTau[1,2] <- 0;\n   tauV[1]<-pow(sigmaV[1],-2);\n",sep="")
  if(modelStructure$Vessel.positiveTows == "zero") Vpos.string = "   for(i in 1:nV) {\n      Vdev[i] <- 0;\n   }\n   vesselTau[1,1] <- 0;\n   vesselTau[1,2] <- 0;\n   sigmaV[1]<-0;\n   tauV[1]<-0;\n"
  
  Vzero.string = ""
  if(modelStructure$Vessel.zeroTows == "fixed") Vzero.string = paste("   pVdev[1] <- 0;\n   for(i in 2:nV) {\n      pVdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselTau[2,1] <- 0;\n","   vesselTau[2,2] <- 0;\n   sigmaV[2]<-0;\n   tauV[2]<-0;\n",sep="")
  if(modelStructure$Vessel.zeroTows == "random") Vzero.string = paste("   for(i in 1:nV) {\n      pVdev[i] ~ dnorm(0,tauV[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselTau[2,1] <- 0;\n","   vesselTau[2,2] <- 0;\n   sigmaV[2] ~ dunif(0,100);\n   tauV[2]<-pow(sigmaV[2],-2);\n",sep="")
  if(modelStructure$Vessel.zeroTows == "random2") Vzero.string = paste("   for(i in 1:nV) {\n      pVdev[i] ~ dnorm(0,tauV[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselTau[2,1] <- 0;\n","   vesselTau[2,2] <- 0;\n   sigmaV[2]<-pow(tauV[2],-1/2);\n   tauV[2]~dgamma(",dgammaNum,",",dgammaNum,");\n",sep="")
  if(modelStructure$Vessel.zeroTows == "randomExpanded") Vzero.string = paste(pVexpanded,"   vesselTau[2,1] <- 0;\n","   vesselTau[2,2] <- 0;\n   tauV[2]<-pow(sigmaV[2],-2);\n",sep="")
  if(modelStructure$Vessel.zeroTows == "zero") Vzero.string = "   for(i in 1:nV) {\n      pVdev[i] <- 0;\n   }\n   vesselTau[2,1] <- 0;\n   vesselTau[2,2] <- 0;\n   sigmaV[2]<-0;\n   tauV[2]<-0;\n"
  
  # combine the strata year interactions into a string
  vessel.string = paste(Vpos.string,Vzero.string)

  if(modelStructure$Vessel.zeroTows == "correlated" & modelStructure$Vessel.positiveTows == "correlated") {
  	# vessel deviations are MVN RE
  	vessel.string = paste("sigmaV[1] <- 0;\n   sigmaV[2] <- 0;\n   vesselTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nV) {\n   vdevs[i,1:2] ~ dmnorm(zs[1:2],vesselTau[1:2,1:2]);\n      Vdev[i] <- min(max(vdevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pVdev[i] <- min(max(vdevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")	
  }
  
  ########################################################################################################################################
  # This section is related to offsets, separate for the positive and binomial models
  # catchability parameter can be set to 1 or estimated as linear or qudratic
  ########################################################################################################################################
  #if(modelStructure$Catchability.positiveTows=="zero") catch.posTows = "   B.pos[1]<-0;\n   B.pos[2] <- 0;\n"
  if(modelStructure$Catchability.positiveTows=="one") catch.posTows = "   logB.pos[1] <- 0;\n   B.pos[1]<-exp(logB.pos[1]);\n   B.pos[2] <- 0;\n"
  if(modelStructure$Catchability.positiveTows=="linear") catch.posTows = "   logB.pos[1] ~ dnorm(0,0.1);\n   B.pos[1]<-exp(logB.pos[1]);\n   B.pos[2] <- 0;\n"
  if(modelStructure$Catchability.positiveTows=="quadratic") catch.posTows = "   logB.pos[1] ~ dnorm(0,0.1);\n   B.pos[1]<-logB.pos[1];\n   logB.pos[2] ~ dnorm(0,0.1);\n   B.pos[2]<-logB.pos[2];\n"
  
  if(modelStructure$Catchability.zeroTows=="zero") catch.zeroTows = "   B.zero[1]<-0;\n   B.zero[2] <- 0;\n"
  if(modelStructure$Catchability.zeroTows=="one") catch.zeroTows = "   logB.zero[1] <- 0;\n   B.zero[1] <- exp(logB.zero[1]);\n   B.zero[2] <- 0;\n"
  if(modelStructure$Catchability.zeroTows=="linear") catch.zeroTows = "   logB.zero[1] ~ dnorm(0,0.1);\n   B.zero[1] <- exp(logB.zero[1]);\n   B.zero[2] <- 0;\n"
  if(modelStructure$Catchability.zeroTows=="quadratic") catch.zeroTows = "   logB.zero[1] ~ dnorm(0,0.1);\n   B.zero[1]<-logB.zero[1];\n   logB.zero[2] ~ dnorm(0,0.1);\n   B.zero[2]<-logB.zero[2];\n"
  
  ########################################################################################################################################
  # This section is related to likelihoods
  # likelihood: lognormal, gamma, inverse gaussian, gamma ECE, lognormal ECE, poisson, negative binomial
  ########################################################################################################################################  
  if(likelihood == "lognormal" | likelihood == "lognormalFixedCV") {
  	# CV is sqrt(exp(sig2)-1) which is approx ~ sigma for sigma is small (< 0.2)
  	likelihood.string = paste("      u.nz[i] <- Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]];\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[1]);\n",sep="")
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]];\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[1]);\n",sep="")
  	}
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   sigma[1] <- sqrt(log(pow(CV[1],2)+1));\n   tau[1] <- pow(sigma[1],-2);\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n",sep="")
  	if(likelihood == "lognormalFixedCV") {
  		prior.string = "   oneOverCV2[1] <- 1;\n   CV[1] <- 1;\n   CV[2] <- 0;\n   sigma[1] <- 1;\n   tau[1] <- 1;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"	
  	}
  }
  
  #2. GAMMA  ##################################################
  if(likelihood == "gamma" | likelihood == "gammaFixedCV") {
  	# gamma in this instance is parameterized in terms of the rate and shape, with mean = a/b, var = a/b2, and CV = 1/sqrt(a)
  	# So parameter 'a' has to be a constant, and 'b' varies by tows - b/c we calculate b = a/mean
    likelihood.string = paste("      u.nz[i] <- exp(min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[1]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[1],b[i]);\n")
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- exp(min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[1]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[1],b[i]);\n")
  	}    
    # for lognormal, gamma prior on 1/sigma2 = 1/CV2. To keep things consistent, gamma prior on a = 1/cv2, b/c CV = 1/sqrt(a)
    # then gamma.b[i] = gamma.a / u[i]
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n",sep="")
  	if(likelihood == "gammaFixedCV") {
  	prior.string = "   oneOverCV2[1] <- 1;\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"  		
  	}	
  }
  
  #3. INVGAUSSIAN  ##################################################  
  if(likelihood == "invGaussian" | likelihood == "invGaussianFixedCV") {
  	# again, parameterize in terms of CV
  	likelihood.string = paste("      u.nz[i] <- exp(min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      lambda[i] <- u.nz[i]*oneOverCV2[1];\n","      scaledLogLike[i] <- -(0.5*log(lambda[i]) - 0.5*logy3[nonZeros[i]] - lambda[i]*pow((y[nonZeros[i]]-u.nz[i]),2)/(2*u.nz[i]*u.nz[i]*y[nonZeros[i]])) + 10000;\n","      zeros.vec[i] ~ dpois(scaledLogLike[i]);\n")
  	
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- exp(min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      lambda[i] <- u.nz[i]*oneOverCV2[1];\n","      scaledLogLike[i] <- -(0.5*log(lambda[i]) - 0.5*logy3[nonZeros[i]] - lambda[i]*pow((y[nonZeros[i]]-u.nz[i]),2)/(2*u.nz[i]*u.nz[i]*y[nonZeros[i]])) + 10000;\n","      zeros.vec[i] ~ dpois(scaledLogLike[i]);\n")
  	} 	
  	
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n",sep="")	
  	if(likelihood == "invGaussianFixedCV") prior.string = "   oneOverCV2[1] <- 1;\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"
  }
 
  #4. POISSON  ##################################################
  if(likelihood == "poisson") {
  	likelihood.string = paste("      u.nz[i] <- exp(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n", "      y[nonZeros[i]] ~ dpois(u.nz[i]);\n",sep="")
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- exp(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n", "      y[nonZeros[i]] ~ dpois(u.nz[i]);\n",sep="")
  	}
  	prior.string = "   oneOverCV2[1] <- 0;\n   CV[1] <- 0;\n   CV[2] <- 0;\n   sigma[1] <- 0;\n   tau[1] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n" 
  }

  #5. NEGATIVE BINOMIAL  ##################################################
  if(likelihood == "negbin") {
  	# calculate p = r / (mean + r)
  	likelihood.string = paste("      u.nz[i] <- r / (r + exp(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]));\n", "      y[nonZeros[i]] ~ dnegbin(u.nz[i],r);\n",sep="")
  	
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- r/(r + exp(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]));\n", "      y[nonZeros[i]] ~ dnegbin(u.nz[i],r);\n",sep="")
  	}
    prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   r <- 1/oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;",sep="")
  }
    
  #6. lognormalECE  ##################################################    
  if(likelihood == "lognormalECE") {
  	# CV is sqrt(exp(sig2)-1) which is approx ~ sigma. Each tow is treated as discrete group (normal, ECE)
  	# if 1, don't do anything 
  	# if 2, multiply by ratio 
  	likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- (G[i]-1)*logratio + (Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[G[i]]);\n",sep="")
  	
  	if(covariates$positive==TRUE) {       
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- (G[i]-1)*logratio + (inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[G[i]]);\n",sep="")
  	} 		
  	
  	# for the ECE model, the normal and extreme distributions each get a separate variance
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   sigma[1] <- sqrt(log(pow(CV[1],2)+1));\n   tau[1] <- pow(sigma[1],-2);\n   oneOverCV2[2] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   CV[2] <- 1/sqrt(oneOverCV2[2]);\n   sigma[2] <- sqrt(log(pow(CV[2],2)+1));\n   tau[2] <- pow(sigma[2],-2);\n   logratio ~ dunif(0,5);\n   ratio <- exp(logratio);\n   alpha.ece[1] <- 1;\n   alpha.ece[2] <- 1;\n   p.ece[1:2] ~ ddirch(alpha.ece[1:2]);\n",sep="")
  }
  
  #7. lognormalECE2  ##################################################    
  if(likelihood == "lognormalECE2") {
  	# CV is sqrt(exp(sig2)-1) which is approx ~ sigma. Each tow is treated as discrete group (normal, ECE)
  	# if 1, don't do anything 
  	# if 2, multiply by ratio 
  	# EW: this implements the Poisson "zeros" trick
  	likelihood.string = paste("      u.nz[i] <- (Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n   scaledLogLike[i] <- -log(p.ece[1]*exp(-logy[nonZeros[i]]-log2pi-log(sigma[1])-pow(logy[nonZeros[i]]-u.nz[i],2)/(2*pow(sigma[1],2))) + p.ece[2]*exp(-logy[nonZeros[i]]-log2pi-log(sigma[2])-pow(logy[nonZeros[i]]-(u.nz[i] + logratio),2)/(2*pow(sigma[2],2)))) + 10000;\n", "      zeros.vec[i] ~ dpois(scaledLogLike[i]);\n",sep="")
  	
  	if(covariates$positive==TRUE) {       
    	# modify the likelihood string to include covariates for the positive tows
  	    likelihood.string = paste("      u.nz[i] <- ((inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n   scaledLogLike[i] <- -log(p.ece[1]*exp(-logy[i]-log2pi-log(sigma[1])-pow(logy[i]-u.nz[i],2)/(2*pow(sigma[1],2))) + p.ece[2]*exp(-logy[i]-log2pi-log(sigma[2])-pow(logy[i]-(u.nz[i] + logratio),2)/(2*pow(sigma[2],2)))) + 10000;\n", "      zeros.vec[i] ~ dpois(scaledLogLike[i]);\n",sep="")    	
  	} 		
  	
  	# for the ECE model, the normal and extreme distributions each get a separate variance
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   sigma[1] <- sqrt(log(pow(CV[1],2)+1));\n   tau[1] <- pow(sigma[1],-2);\n   oneOverCV2[2] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   CV[2] <- 1/sqrt(oneOverCV2[2]);\n   sigma[2] <- sqrt(log(pow(CV[2],2)+1));\n   tau[2] <- pow(sigma[2],-2);\n   logratio ~ dunif(0,5);\n   ratio <- exp(logratio);\n   alpha.ece[1] <- 1;\n   alpha.ece[2] <- 1;\n   p.ece[1:2] ~ ddirch(alpha.ece[1:2]);\n",sep="")
  }
  
  #8. gammaECE  ##################################################
  if(likelihood == "gammaECE") {
  	# gamma in this instance is parameterized in terms of the rate and shape, with mean = a/b, var = a/b2, and CV = 1/sqrt(a)
  	# So parameter 'a' has to be a constant, and 'b' varies by tows - b/c we calculate b = a/mean
    likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- exp((G[i]-1)*logratio + min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[G[i]]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[G[i]],b[i]);\n")
      
  	if(covariates$positive==TRUE) {
  	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- exp((G[i]-1)*logratio + min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[G[i]]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[G[i]],b[i]);\n")
  	} 		    
      
      # for lognormal, gamma prior on 1/sigma2 = 1/CV2. To keep things consistent, gamma prior on a = 1/cv2, b/c CV = 1/sqrt(a)
      # then gamma.b[i] = gamma.a / u[i]
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   oneOverCV2[2] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   gamma.a[2] <- oneOverCV2[2];\n   CV[2] <- 1/sqrt(oneOverCV2[2]);\n   logratio ~ dunif(0,5);\n   ratio <- exp(logratio);\n   alpha.ece[1] <- 1;\n   alpha.ece[2] <- 1;\n   p.ece[1:2] ~ ddirch(alpha.ece[1:2]);\n",sep="")
  }
      
  #9. gammaECE  ##################################################
  if(likelihood == "gammaECE2") {
  	# gamma in this instance is parameterized in terms of the rate and shape, with mean = a/b, var = a/b2, and CV = 1/sqrt(a)
  	# So parameter 'a' has to be a constant, and 'b' varies by tows - b/c we calculate b = a/mean
    likelihood.string = paste("      u.nz[i] <- exp(min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      scaledLogLike[i] <- -log(p.ece[1]*exp(gamma.a[1]*log(gamma.a[1]/u.nz[i]) + (gamma.a[1]-1)*logy[nonZeros[i]] -(gamma.a[1]/u.nz[i])*y[nonZeros[i]] - loggam(gamma.a[1])) + p.ece[2]*exp(gamma.a[2]*log(gamma.a[2]/(u.nz[i]*exp(logratio))) + (gamma.a[2]-1)*logy[nonZeros[i]] -(gamma.a[2]/(u.nz[i]*exp(logratio)))*y[nonZeros[i]] - loggam(gamma.a[2]))) + 10000;\n",     "zeros.vec[i] ~ dpois(scaledLogLike[i]);\n",sep="")
      
  	if(covariates$positive==TRUE) {
  	# modify the likelihood string to include covariates for the positive tow
  	likelihood.string = paste("      u.nz[i] <- exp(min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + Vdev[vessel[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      scaledLogLike[i] <- -log(p.ece[1]*exp(gamma.a[1]*log(gamma.a[1]/u.nz[i]) + (gamma.a[1]-1)*logy[i] -(gamma.a[1]/u.nz[i])*y[i] - loggam(gamma.a[1])) + p.ece[2]*exp(gamma.a[2]*log(gamma.a[2]/(u.nz[i]*exp(logratio))) + (gamma.a[2]-1)*logy[i] -(gamma.a[2]/(u.nz[i]*exp(logratio)))*y[i] - loggam(gamma.a[2]))) + 10000;\n",     "zeros.vec[i] ~ dpois(scaledLogLike[i]);\n")   	 	
  	}
    # for lognormal, gamma prior on 1/sigma2 = 1/CV2. To keep things consistent, gamma prior on a = 1/cv2, b/c CV = 1/sqrt(a)
    # then gamma.b[i] = gamma.a / u[i]
  	prior.string = paste("   oneOverCV2[1] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   oneOverCV2[2] ~ dgamma(",dgammaNum,",",dgammaNum,");\n   gamma.a[2] <- oneOverCV2[2];\n   CV[2] <- 1/sqrt(oneOverCV2[2]);\n   logratio ~ dunif(0,5);\n   ratio <- exp(logratio);\n   alpha.ece[1] <- 1;\n   alpha.ece[2] <- 1;\n   p.ece[1:2] ~ ddirch(alpha.ece[1:2]);\n",sep="")
  }
  
  ####################################################################
  # This section is related to year deviations, default is to make them uncorrelated
  # if strata-year effects are estimated as fixed effects, parameters are redundant, so year deviations set to 0
  str1 = paste("      Ydev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "fixed") {str1 = "      Ydev[i] <- 0;\n"}
  str2 = paste("      pYdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "fixed") {str2 = "      pYdev[i] <- 0;\n"}
  
  year.dev.string = paste("   for(i in 1:nY) { # year deviations, always fixed \n",str1,str2,"}\n   yearTau[1,1] <- 0;\n   yearTau[1,2] <- 0;\n   yearTau[2,1]<-0;\n   yearTau[2,2] <- 0;\n",sep="")
  if(modelStructure$year.deviations == "correlated") {
  	year.dev.string = paste("   yearTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nY) {\n   ydevs[i,1:2] ~ dmnorm(zs[1:2],yearTau[1:2,1:2]);\n      Ydev[i] <- min(max(ydevs[i,1],",logBounds[1],"),",logBounds[1],");\n      pYdev[i] <- min(max(ydevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")
  }
  
  ####################################################################
  # This section is related to strata deviations, default is to make them uncorrelated
  # if strata-year effects are estimated as fixed effects, parameters are redundant, so strata deviations set to 0
  str1 = paste("      Sdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "fixed") {str1 = "      Sdev[i] <- 0;\n"}
  str2 = paste("      pSdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "fixed") {str2 = "      pSdev[i] <- 0;\n"}
  strata.dev.string = paste("   for(i in 2:nS) { # strata deviations, always fixed \n", str1, str2,"}\n   strataTau[1,1] <- 0;\n   strataTau[1,2] <- 0;\n   strataTau[2,1]<-0;\n   strataTau[2,2] <- 0;\n")
  if(modelStructure$strata.deviations == "correlated") {
  	strata.dev.string = paste("   strataTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 2:nS) {\n   sdevs[i,1:2] ~ dmnorm(zs[1:2],strataTau[1:2,1:2]);\n      Sdev[i] <- min(max(sdevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pSdev[i] <- min(max(sdevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")
  }
  
  # This set of if() blocks is used to include optional strings for covariates...these are all the priors, vague normal
  covarString = ""
  if(covariates$binomial == FALSE) {
  	covarString = paste(covarString, "	C.bin[1] <- 0;\n",sep="")	
  } 
  if(covariates$binomial == TRUE) {
  	covarString = paste(covarString, "	for(i in 1:nX.binomial) {\n","		C.bin[i] ~ dnorm(0,0.001);\n","	}\n",sep="")		
  }
  if(covariates$positive == FALSE) {
  	covarString = paste(covarString, "	C.pos[1] <- 0;\n",sep="")	
  } 
  if(covariates$positive == TRUE) {
  	covarString = paste(covarString, "	for(i in 1:nX.pos) {\n","		C.pos[i] ~ dnorm(0,0.001);\n","	}\n",sep="")		
  }
  
  # This IF statement is just to modify the expected value of the binomial model to include covariates IF they are specified 
  if(pres_link == "logit") logit.p.string = "logit(p.z[i]) <- pVdev[vessel[i]] + pVYdev[vesselYear[i]] + pSdev[strata[i]] + pYdev[year[i]] + pSYdev[strataYear[i]] + B.zero[1]*effort[i] + B.zero[2]*effort2[i];\n"
  if(pres_link == "probit") logit.p.string = "probit(p.z[i]) <- pVdev[vessel[i]] + pVYdev[vesselYear[i]] + pSdev[strata[i]] + pYdev[year[i]] + pSYdev[strataYear[i]] + B.zero[1]*effort[i] + B.zero[2]*effort2[i];\n"
  if(covariates$binomial == TRUE) {
  	if(pres_link == "logit") logit.p.string = "logit(p.z[i]) <- inprod(C.bin[1:nX.binomial],X.bin[i,]) + pVdev[vessel[i]] + pVYdev[vesselYear[i]] + pSdev[strata[i]] + pYdev[year[i]] + pSYdev[strataYear[i]] + B.zero[1]*effort[i] + B.zero[2]*effort2[i];\n"	
  	if(pres_link == "probit") logit.p.string = "probit(p.z[i]) <- inprod(C.bin[1:nX.binomial],X.bin[i,]) + pVdev[vessel[i]] + pVYdev[vesselYear[i]] + pSdev[strata[i]] + pYdev[year[i]] + pSYdev[strataYear[i]] + B.zero[1]*effort[i] + B.zero[2]*effort2[i];\n"	
  }
  
  deltaGLM = 
  paste(
  "
  model {\n",
     prior.string,
     covarString,	
     
  "   # next deal with catchability parameters\n",
     catch.posTows,
     catch.zeroTows,
  "  # these zeros are optionally for MVN random effects
     zs[1] <- 0; 
     zs[2] <- 0;
     # first stratum is set to 0 for identifiability
     Sdev[1] <- 0;
     pSdev[1] <- 0;\n
  ",
  year.dev.string,
  strata.dev.string,   
  stratayear.string,   
  vesselyear.string,
  vessel.string,   	    	    
  "  # for each group, calculate the probability of zero tows and the mean
     # evaluate the likelihood  of non-zero trawls   
     for(i in 1:nNonZeros) {\n",
     	likelihood.string,
  "   }
     
     # evaluate the likelihood  of zero / non-zero trawls
     for(i in 1:n) {
     	 # group represents which vessel x strata combo\n",
  logit.p.string,
  "     isNonZeroTrawl[i] ~ dbern(p.z[i]);
     }
  }
  ",sep="")
  # write this to text file
  if(write.model) cat(deltaGLM, file = model.name)
  
  # fit the model and return the model object
  modelFit = NA
  
  if(fit.model) {
jags.params=c("Ydev","Sdev","SYdev","VYdev","pYdev","pSdev","pSYdev","pVYdev","B.zero","B.pos","sigmaSY","sigmaVY","sigmaV","CV","ratio","p.ece","yearTau","strataTau","strataYearTau","vesselYearTau","C.pos","C.bin","Vdev","pVdev")
    jags.data = list("y","logy3","logy","log2pi","effort","effort2","logeffort", "logeffort2","nonZeros", "n", "isNonZeroTrawl", "nNonZeros", "nV","nVY", "nSY", "nS", "nY", "year", "vesselYear", "vessel","strata", "strataYear","zeros.vec","R","X.pos","X.bin","nX.binomial","nX.pos")
  
    capture.output(jags.data, file="jags_data.txt")
    if(Parallel==TRUE) {
    	modelFit= jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.name, n.chains = mcmc.control$chains, n.burnin = mcmc.control$burn, n.thin = mcmc.control$thin, n.iter = as.numeric(mcmc.control$burn+mcmc.control$iterToSave), DIC = TRUE)
    }
    if(Parallel==FALSE) {
  	  modelFit = jags(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.name, n.chains = mcmc.control$chains, n.burnin = mcmc.control$burn, n.thin = mcmc.control$thin, n.iter = as.numeric(mcmc.control$burn+mcmc.control$iterToSave), DIC = TRUE)
    }
  }
  
  ######################################################################
  # Create a list object that contains all of the properties of this run
  ######################################################################
  # modelStructure (list)
  # model.name (string)
  # fit.model (boolean)
  # mcmc.control (list)
  # Parallel (boolean)
  ######################################################################
  functionCall = list("modelStructure"=modelStructure, "model.name"=model.name, "fit.model"=fit.model, "mcmc.control"=mcmc.control, "Parallel"=Parallel,"likelihood"=likelihood,"covariates"=covariates) # species, likelihood, parameters TRUE/FALSE
  
  ######################################################################
  # Return a list of which parameters are actually estimated
  ######################################################################
  estimatedParameters = NA
  if(fit.model) {
    estimatedParameters = data.frame("Parameter" = colnames(modelFit$BUGSoutput$sims.matrix))
    estimatedParameters$Estimated = rep(TRUE, dim(modelFit$BUGSoutput$sims.matrix)[2])
    vars = apply(modelFit$BUGSoutput$sims.matrix,2,var) # figure out which cols have var = 0
    estimatedParameters$Estimated[which(vars==0)] = FALSE
  }
  
  return(c(modelFit,functionCall,estimatedParameters,Species=Species,Data=Data))

}
