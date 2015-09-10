########################################################
####### This block of code is related to processing output
doMCMCDiags = function(datalist, strata.limits=strata.limits, directory, mods, StrataWeights="StrataAreas", McmcDiagnostics=FALSE) {

  # Load tagged list of data
  attach(datalist)
  on.exit( detach(datalist) )  
  
  # Load data locally
  attach(Data)
  on.exit( detach(Data), add = TRUE )

  # Identify strata and year for StratYear values
  StrataTable = data.frame(strataYear = levels(strataYear), 
        strata = sapply(levels(strataYear), FUN = function(Char) {
                   Split <- strsplit(Char, ":")[[1]]; paste(Split[-length(Split)], collapse=":")}), 
          year = sapply(levels(strataYear), FUN = function(Char) {
                    Split <- strsplit(Char, ":")[[1]]; Split[length(Split)]}), 
                    Area_Hectares = rep(NA, nlevels(strataYear)))
  for(i in 1:nrow(StrataTable)){
    Row = which(strata.limits[,'STRATA']==StrataTable[i,'strata'])
    StrataTable[i,'Area_Hectares'] = sum(SA3[SA3[,'MAX_LAT_DD']<=strata.limits[Row,'NLat'] & SA3[,'MIN_LAT_DD']>=strata.limits[Row,'SLat'] & SA3[,'MIN_DEPTH_M']>=strata.limits[Row,'MinDepth'] & SA3[,'MAX_DEPTH_M']<=strata.limits[Row,'MaxDepth'],'AREA_HECTARES'])
  }

  # Make folder for plots
  Species = species = mods[[1]]$Species
  SpeciesFolder = paste(getwd(),"/",Species,"_FinalDiagnostics/",sep="")
  dir.create(SpeciesFolder, showWarnings=FALSE)

  # Make objects for saving output
  Indices = array(NA, dim=c(nlevels(year),length(mods),2))
  IndicesByStrata = array(NA, dim=c(nlevels(year),nlevels(strata),length(mods),2))

  ######################
  # Display data
  ######################

  # Plot data by year, depth, and latitude
  PlotData(Data=Data, FileName="", Folder=SpeciesFolder)
  # Plot location of data
  MapData(Data=Data, strata.limits=strata.limits, SA=SA3, Folder=SpeciesFolder)
  # Save mods for later usage
  Save = list(mods=mods, Data=Data)
  save(Save, file=paste(SpeciesFolder,"Save.RData",sep=""))
  # Load old mods (if you want to)
  #load(file=paste(SpeciesFolder,"Save.RData",sep=""))
  #attach(Save)

  # Format data
  out <- list()
  for(ModelNumber in 1:length(mods)){

    # Make folder
    Folder = paste(SpeciesFolder,"/Model=",ModelNumber,"/",sep="")
    dir.create(Folder, showWarnings=FALSE)

    # Unpack results
    Model = mods[[ModelNumber]]
    McmcArray = as.array(Model$BUGSoutput$sims.array)
    McmcList = vector("list",length=dim(McmcArray)[2])
    for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
    McmcList = mcmc.list(McmcList)
    BugsList = Model$BUGSoutput$sims.list

    # Record details
    capture.output(Model$mcmc.control, file=paste(Folder,"mcmc_control.txt",sep=""))
    capture.output(list(Model$likelihood,Model$modelStructure), file=paste(Folder,"Model_Structure.txt",sep=""))
    capture.output(Model$BUGSoutput$summary, file=paste(Folder,"BUGSoutput_summary.txt",sep=""))
    capture.output(Model$model, file=paste(Folder,"deltaGLM.txt",sep=""))
    save(Model, file=paste(Folder,"mods_for_this_run.RData",sep=""))

    # Load parToMinitor from Eric's code
    parToMonitor = which(Model$Estimated)
    parnames = Model$Parameter[parToMonitor]

    ######################
    # Convergence diagnostics
    ######################

    # TracePlot
    # maxDims=10; Nkeep=5000; parToMonitor=parToMonitor; parnames=parnames[parToMonitor]; FileName=""; Type="Trace"; Folder=Folder
    ConvergencePlot(McmcArray, maxDims=8, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="Trace", Folder=Folder, Model=Model)
    # ACF plot
    ConvergencePlot(McmcArray, maxDims=8, Nkeep=500, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="ACF", Folder=Folder, Model=Model)
    # Variance density plots
    if(length(grep("sigma",parnames)>0)){
      ConvergencePlot(McmcArray, maxDims=8, Nkeep=500, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="VarianceDensity", Folder=Folder, Model=Model)
    }
    # Correlation density plots
    #if(Model$modelStructure$StrataYear.positiveTows == "correlated" || Model$modelStructure$VesselYear.positiveTows == "correlated") ConvergencePlot(McmcArray, maxDims=8, Nkeep=500, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="CorrelationDensity", Folder=Folder)

    # Convergence statistics -- This doesn't always converge for either Geweke or Gelman-Rubin statistics
    if(McmcDiagnostics==TRUE){
      McmcDiagnosticsPlot(McmcList, parToMonitor, FileName="", Folder=Folder, Geweke=FALSE)
    }

    #######################
    # Model fit diagnostics
    #######################

    # Visualize the realized offset
    PlotOffset(Data=Data, BugsList=BugsList, FileName="", Folder=Folder)
    # Posterior predictive distribution for positive catches
    WAIC <- PosteriorPredictive(Data=Data, Model=Model, FileName="", Folder=Folder)
    # JAGS indices of abundance
    McmcIndices = ComputeIndices(Data=Data, Model=Model, FileName="", Folder=Folder, Weights=StrataWeights, StrataTable=StrataTable)
    # MLE indices of abundance
    MleIndices = try(ComputeMleIndices(Data=Data, Model=Model, FileName="", Folder=Folder, Weights=StrataWeights, StrataTable=StrataTable, Run=TRUE), silent=TRUE)
    if(inherits(MleIndices, "try-error")==TRUE){
      MleIndices = ComputeMleIndices(Data=Data, Model=Model, FileName="", Folder=Folder, Weights=StrataWeights, StrataTable=StrataTable, Run=FALSE)
    }
    out[[ifelse(is.null(names(mods)), ModelNumber, names(mods)[ModelNumber])]] <- list(MCMC.byYear = McmcIndices$byYear, 
      MCMC.byYearAndStrata = McmcIndices$byYearAndStrata, MLE.byYear = MleIndices$byYear, MLE.byYearAndStrata = MleIndices$byYearAndStrata, WAIC = WAIC)

    # Compare JAGS and MLE
    jpeg(paste(Folder,"/","","Index_Comparison.jpg",sep=""),width=2*3,height=2*3,res=200,units="in")
    par(mfrow=c(2,2), mar=c(2.5,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    matplot(cbind(McmcIndices$byYearAndStrata[,c('PresMedian','RawPres')], MleIndices$byYearAndStrata$Pres), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="StrataYear Presence", ylim=c(0,max(cbind(McmcIndices$byYearAndStrata[,c('PresMedian','RawPres')], MleIndices$byYearAndStrata$Pres),na.rm=TRUE)))
    matplot(cbind(McmcIndices$byYearAndStrata[,c('PosMedian','RawPos')], MleIndices$byYearAndStrata$Pos), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="StrataYear Positive catch", ylim=c(0,max(cbind(McmcIndices$byYearAndStrata[,c('PosMedian','RawPos')], MleIndices$byYearAndStrata$Pos),na.rm=TRUE)))
    matplot(cbind(McmcIndices$byYearAndStrata[,c('IndexMedian','Raw')], MleIndices$byYearAndStrata$Index), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="StrataYear Index", ylim=c(0,max(cbind(McmcIndices$byYearAndStrata[,c('IndexMedian','Raw')], MleIndices$byYearAndStrata$Index),na.rm=TRUE)))
    matplot(cbind(McmcIndices$byYear[,c('IndexMedian','Raw')], MleIndices$byYear$Index), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="Year Index", ylim=c(0,max(cbind(McmcIndices$byYear[,c('IndexMedian','Raw')], MleIndices$byYear$Index),na.rm=TRUE)))
    dev.off()

    # Make easier to read version
    jpeg(paste(Folder,"/","","Index_with_95CI.jpg",sep=""),width=4,height=4,res=200,units="in")
    par(mfrow=c(1,1), mar=c(2.5,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    plot(McmcIndices$byYear[,'IndexMedian'], col="black", xlab="Year", ylab="Index",main="Year Index", ylim=c(0,max(McmcIndices$byYear[,'IndexMedian']*exp(1.96*McmcIndices$byYear[,'SdLog']),na.rm=TRUE)))
    for(i in 1:nrow(McmcIndices$byYear)) lines( x=rep(i,2), y=McmcIndices$byYear[i,'IndexMedian']*exp(c(-1.96,1.96)*McmcIndices$byYear[i,'SdLog']), col="red")
    dev.off()

    # Save McmcIndices and CV
    Indices[,ModelNumber,1] = McmcIndices$byYear$IndexMean
    Indices[,ModelNumber,2] = McmcIndices$byYear$SdLog
    for(StratI in 1:nlevels(strata)){
      Which = which(McmcIndices$byYearAndStrata$Strata==levels(McmcIndices$byYearAndStrata$Strata)[StratI])
      IndicesByStrata[,StratI,ModelNumber,1] = McmcIndices$byYearAndStrata$IndexMean[Which]
      IndicesByStrata[,StratI,ModelNumber,2] = McmcIndices$byYearAndStrata$SdLog[Which]
    }
  }

  # Plot Index and CV for all model configurations
  jpeg(paste(SpeciesFolder,"Index_Comparison.jpg",sep=""),width=1*4,height=2*4,res=200,units="in")
  par(mfrow=c(2,1), mgp=c(1.25,0.25,0), mar=c(3,3,1,0), tck=-0.02)
  matplot(Indices[,,1], col="black", lty="solid", type="b", xlab="Year", ylab="Biomass index", ylim=c(0,max(Indices[,,1],na.rm=T)))
  matplot(Indices[,,2], col="black", lty="solid", type="b", xlab="Year", ylab="Index CV", ylim=c(0,max(Indices[,,2],na.rm=T)))
  dev.off()

  # Plot Index and CV | Strata for all model configurations
  Log = ""
  jpeg(paste(SpeciesFolder,"Index_Comparison_by_strata.jpg",sep=""),width=2*3,height=nlevels(strata)*3,res=200,units="in")
  par(mfrow=c(nlevels(strata),2), mgp=c(1.25,0.25,0), mar=c(3,3,1,0), tck=-0.02, oma=c(0,2,2,0))
  for(StratI in 1:nlevels(strata)){
    matplot(IndicesByStrata[,StratI,,1], col="black", log=Log, lty="solid", type="b", xlab="Year", ylab="Biomass index", ylim=list( c(0,max(IndicesByStrata[,,,1])), range(IndicesByStrata[,,,1]) )[[ifelse(Log=="",1,2)]])
    if(StratI==1) mtext(side=3, outer=FALSE, line=1, text="Biomass index", cex=1.5)
    mtext(side=2, outer=FALSE, line=2, text=levels(strata)[StratI], cex=1.5)
    matplot(IndicesByStrata[,StratI,,2], col="black", lty="solid", type="b", xlab="Year", ylab="Index CV", ylim=c(0,max(IndicesByStrata[,,,2],na.rm=T)))
    if(StratI==1) mtext(side=3, outer=FALSE, line=1, text="Index CV", cex=1.5)
  }
  dev.off()
  invisible(out)
}
