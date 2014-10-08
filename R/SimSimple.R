
SimSimple = function(MeanEncounter, Nyears=10, Nstrata=15, Nvessels=4, Obsperyear=175, sigmaV=rep(1,2), sigmaVY=rep(1,2), sigmaS=rep(1,2), sigmaY=rep(1,2), sigmaSY=rep(1,2), sigmaResid=1){

  # Simulate effects
  betaY = array( rnorm(Nyears*2,sd=sigmaY), dim=c(2,Nyears))
  betaS = array( rnorm(Nstrata*2,sd=sigmaY), dim=c(2,Nstrata))
  betaSY = array( rnorm(Nstrata*Nyears*2,sd=sigmaSY), dim=c(2,Nstrata,Nyears))
  betaV = array( rnorm(Nvessels*2,sd=sigmaV), dim=c(2,Nvessels))
  betaVY = array( rnorm(Nvessels*Nyears*2,sd=sigmaVY), dim=c(2,Nvessels,Nyears))

  # Simulate data
  DF = data.frame()
  for(YearI in 1:Nyears){
    Data = cbind( "YEAR"=rep(YearI,Obsperyear), "STRATUM"=sample(1:Nstrata,replace=TRUE,size=Obsperyear), "VESSEL"=sample(1:Nvessels,replace=TRUE,size=Obsperyear), "AREA_SWEPT_HA"=rep(1,Obsperyear) )
    Data = cbind( Data, "SPECIES_WT_KG"=rbinom( n=Obsperyear, size=1, prob=plogis(qlogis(MeanEncounter) + betaY[1,][Data[,'YEAR']] + betaS[1,][Data[,'STRATUM']] + betaSY[1,,][Data[,c('STRATUM','YEAR')]] + betaV[1,][Data[,'VESSEL']] + betaVY[1,,][Data[,c('VESSEL','YEAR')]]) ) )
    alpha = 1 / sigmaResid^2
    beta = 1 / ( sigmaResid^2 * exp(betaY[2,][Data[,'YEAR']] + betaS[2,][Data[,'STRATUM']] + betaSY[2,,][Data[,c('STRATUM','YEAR')]] + betaV[2,][Data[,'VESSEL']] + betaVY[2,,][Data[,c('VESSEL','YEAR')]]) )
    Data[,'SPECIES_WT_KG'] = ifelse( Data[,'SPECIES_WT_KG']==1, rgamma(Obsperyear, shape=alpha, rate=beta), 0 )
    DF = rbind(DF, Data)
  }
  DF[,'YEAR'] = DF[,'YEAR'] + 2000
  # Re-organize stratum details
  DF = as.data.frame(DF)
  names(DF) = colnames(Data)
  DF = data.frame(DF, "BEST_DEPTH_M"=rep(1,Obsperyear), "BEST_LAT_DD"=DF[,'STRATUM'] )
  
  Return = list( "DF"=DF, "betaY"=betaY, "betaS"=betaS, "betaSY"=betaSY, "betaV"=betaV, "betaVY"=betaVY )
  return(Return)
}
