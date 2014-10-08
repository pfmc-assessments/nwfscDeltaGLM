ComputeIndices = function(Data, Model, FileName, maxDims=6, Folder=NA, Weights="StrataAreas", StrataTable, PlotStrataYearMcmc=TRUE){
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff -- listed by search()
  attach(Model$BUGSoutput$sims.list)
  #attach(Data)
  modelStructure = Model$modelStructure
  Dist = Model$likelihood
  
  # Compute average AreaSwept - This AreaSwept used to calculate densities, and matters if Catch is a nonlinear function of AreaSwept
  MeanLogAreaSwept = mean(Data$logeffort)
  
  # New objects
  Chains = array(NA,dim=c(nrow(Sdev),length(unique(StrataTable[,'year'])),length(unique(StrataTable[,'strata'])),2))
  Year = Strata = Area = PosMedian = PresMedian = IndexMedian = IndexMedianWeighted = PosMean = PresMean = IndexMean = IndexMeanWeighted = CvMedian = SdLog = RawPos = RawPres = Raw = RawWeighted = RawVar = RawVarWeighted = RawCV = matrix(NA,nrow=length(unique(StrataTable[,'year'])),ncol=length(unique(StrataTable[,'strata']))) 
  CvMedianYear = SdLogYear = rep(NA,length(unique(StrataTable[,'year'])))
  
  # Calculate values
  for(YearI in 1:nrow(PosMean)){
    for(StratI in 1:ncol(PosMean)){
      # Derived indicators
      StrataYearI = which(levels(strataYear)==paste(toupper(levels(strata)[StratI]),":",levels(year)[YearI],sep=""))
      AreaI = which(StrataTable[,'strataYear']==paste(levels(strata)[StratI],":",levels(year)[YearI],sep=""))
      Which = which(strata==levels(strata)[StratI] & year==levels(year)[YearI])
      # Year, strata, and area
      Year[YearI,StratI] = levels(year)[YearI]
      Strata[YearI,StratI] = levels(strata)[StratI]
      if(Weights=="StrataAreas") Area[YearI,StratI] = StrataTable[AreaI,'Area_Hectares']
      if(Weights=="Equal") Area[YearI,StratI] = 1
      # Save Positive catch chain in normal-space and correct for transformation biases
      Chains[,YearI,StratI,1] = exp( cMx(Sdev)[,StratI] + cMx(Ydev)[,YearI] + cMx(SYdev)[,StrataYearI] + log(1)*B.pos[,1] + log(1)^2*B.pos[,2] ) # wardJAGS uses logeffort offset
        # Lognormal -- Bias correction 
        if(Dist=="lognormal"){
          Sigma = sqrt(log(CV[,1]^2+1))        
          Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1] * exp(Sigma^2/2)
        }
        # LognormalECE -- Bias correction + incorporate ECE
        if(Dist=="lognormalECE"){
          Sigma = sqrt(log(CV[,1]^2+1))        
          Sigma2 = sqrt(log(CV[,2]^2+1))        
          Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1]*exp(Sigma^2/2)*p.ece[,1] + ratio*Chains[,YearI,StratI,1]*exp(Sigma2^2/2)*p.ece[,2]
        }
        # GammaECE -- incorporate ECE 
        if(Dist=="gammaECE"){
          Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1]*p.ece[,1] + ratio*Chains[,YearI,StratI,1]*p.ece[,2]
        }
        # Don't make mean-unbiased for unobserved vessel
        #if(modelStructure$VesselYear.positiveTows=="random") Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1] * exp(sigmaVY[,1]^2/2)
      # Save Presence/Absence chain in normal-space and correct for transformation biases
      Chains[,YearI,StratI,2] = plogis( cMx(pSdev)[,StratI] + cMx(pYdev)[,YearI] + cMx(pSYdev)[,StrataYearI] + (1)*B.zero[,1] + (1)^2*B.zero[,2] )   # wardJAGS predicts the probability of 0 catch, and uses offset as effort, not logeffort
        # Don't make mean-unbiased for unobserved vessel (this was done incorrectly anyway)
        #if(modelStructure$VesselYear.zeroTows=="random") Chains[,YearI,StratI,2] = mean(plogis(rnorm(n=dim(Chains)[1]*1e3, mean=qlogis(Chains[,YearI,StratI,2]), sd=sigmaVY[,2])))
      # Index (median)
      PosMedian[YearI,StratI] = median( Chains[,YearI,StratI,1] )  # / 2e4 # Convert kilograms to metric tons
      PresMedian[YearI,StratI] = median( Chains[,YearI,StratI,2] )   
      IndexMedian[YearI,StratI] = median( Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2] ) 
      IndexMedianWeighted[YearI,StratI] = IndexMedian[YearI,StratI] * Area[YearI,StratI]
      # Index (mean)
      PosMean[YearI,StratI] = mean( Chains[,YearI,StratI,1] )  # / 2e4 # Convert kilograms to metric tons
      PresMean[YearI,StratI] = mean( Chains[,YearI,StratI,2] )   
      IndexMean[YearI,StratI] = mean( Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2] ) 
      IndexMeanWeighted[YearI,StratI] = IndexMean[YearI,StratI] * Area[YearI,StratI]
      # CV of median (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r)
      Temp = Area[YearI,StratI] * Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2]
      CvMedian[YearI,StratI] = sqrt(var(Temp)) / median(Temp)
      # CV if median (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r)
      Temp = Area[YearI,StratI] * Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2]
      SdLog[YearI,StratI] = sd(log(Temp))
      # Raw
      RawPos[YearI,StratI] = mean(ifelse(Data[Which,'HAUL_WT_KG']>0,Data[Which,'HAUL_WT_KG']/Data[Which,'effort'],NA),na.rm=TRUE) 
      RawPres[YearI,StratI] = mean(Data[Which,'HAUL_WT_KG']>0)    
      Raw[YearI,StratI] = RawPos[YearI,StratI] * RawPres[YearI,StratI]
      RawWeighted[YearI,StratI] = Raw[YearI,StratI] * Area[YearI,StratI]
      RawVar[YearI,StratI] = var(Data[Which,'HAUL_WT_KG']/Data[Which,'effort'],na.rm=TRUE) / length(Which) 
      RawVarWeighted[YearI,StratI] = RawVar[YearI,StratI] * Area[YearI,StratI]^2 
      RawCV[YearI,StratI] = sqrt( RawVarWeighted[YearI,StratI] ) / RawWeighted[YearI,StratI] 
    }
    # CV of median (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r")
    Temp = Area[YearI,StratI] * rowSums( cMx(Chains[,YearI,,1]) * cMx(Chains[,YearI,,2]) )
    CvMedianYear[YearI] = sqrt(var(Temp)) / median(Temp)
    # SD of log of index (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r")
    Temp = Area[YearI,StratI] * rowSums( cMx(Chains[,YearI,,1]) * cMx(Chains[,YearI,,2]) )
    SdLogYear[YearI] = sd(log(Temp))
  } # 1085-115

  # Plot MCMC chains for each Strata:Year
  if(PlotStrataYearMcmc == TRUE){
    for(Type in c("Trace","ACF")){
      for(ChainI in 1:2){
        Nstrat = length(unique(strataYear[nonZeros]))
        Nplots = ceiling( Nstrat / maxDims^2 )
        for(PlotI in 1:Nplots){
          ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
          Ncol=ceiling(sqrt(length(ParSet)))
          Nrow = ceiling(length(ParSet)/Ncol)
          jpeg(paste(Folder,"/",FileName,"Chain_",c("Positive","Presence")[ChainI],"_",Type,"_by_StrataYear_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=150,units="in")
            par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
            for(StrataYearI in ParSet){
              StrataI = strata[which(strataYear==unique(strataYear)[StrataYearI])[1]]
              YearI = year[which(strataYear==unique(strataYear)[StrataYearI])[1]]
              Mat = matrix(Chains[,YearI,StrataI,ChainI],ncol=Model$mcmc.control$chains,byrow=FALSE)
              if(Type=="Trace"){
                matplot(Mat, type="l", lty="solid", col=rainbow(Model$mcmc.control$chains,alpha=0.4), main=paste(YearI,StrataI), xaxt="n", xlab="",ylab="")
              }
              if(Type=="ACF"){
                Acf = apply(Mat,MARGIN=2,FUN=function(Vec){acf(Vec,plot=FALSE)$acf})
                matplot(Acf, type="h", lty="solid", col=rainbow(Model$mcmc.control$chains,alpha=0.4), main=paste(YearI,StrataI), xaxt="n", xlab="",ylab="", ylim=c(0,1), lwd=2)
              }
            }
          dev.off()
        }
      }
    }
  }
  
  # Compile into matrices
  Results1 = data.frame(Year=as.vector(Year), Strata=as.vector(Strata), Raw=as.vector(RawWeighted), RawCV=as.vector(RawCV), IndexMedian=as.vector(IndexMedianWeighted), IndexMean=as.vector(IndexMeanWeighted), CvMedian=as.vector(CvMedian), SdLog=as.vector(SdLog), Area=as.vector(Area), PosMedian=as.vector(PosMedian), PresMedian=as.vector(PresMedian), PosMean=as.vector(PosMean), PresMean=as.vector(PresMean), RawPos=as.vector(RawPos), RawPres=as.vector(RawPres), RawSD=ifelse(as.vector(RawVarWeighted)==0,NA,as.vector(sqrt(RawVarWeighted))))
  Results2 = data.frame(Year=Year[,1], Raw=rowSums(RawWeighted,na.rm=TRUE), RawCV=sqrt(rowSums(RawVarWeighted,na.rm=TRUE))/rowSums(RawWeighted,na.rm=TRUE), IndexMedian=rowSums(IndexMedianWeighted,na.rm=TRUE), IndexMean=rowSums(IndexMeanWeighted,na.rm=TRUE), CvMedian=CvMedianYear, SdLog=SdLogYear)
  
  # Detach stuff -- listed by search()
  #detach(Data)
  detach(Model$BUGSoutput$sims.list)
  
  # Write and print output
  write.csv(Results1,file=paste(Folder,"/",FileName,"ResultsByYearAndStrata.csv",sep=""))
  write.csv(Results2,file=paste(Folder,"/",FileName,"ResultsByYear.csv",sep=""))

  # Return output
  Return = list(Results1=Results1, Results2=Results2)
  return(Return)
}

