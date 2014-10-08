############################################################################
# FUNCTION processData WRITTEN BY JIM THORSON & ERIC WARD, UPDATED 9/30/2012. 
# EMAIL: JAMES.THORSON@NOAA.GOV, ERIC.WARD@NOAA.GOV
############################################################################
processData = function(Truncate=0) {
  
  print("Necessary column names for masterDat:")
  print("1. BEST_DEPTH_M -> tow depth in meters")
  print("2. BEST_LAT_DD -> tow latitude in degrees ")
  print("3. [insert species name] -> tow catch in kilograms ")
  print("4. YEAR -> calendar year, or time-strata")
  print("5. AREA_SWEPT_MSQ -> area-swept in square meters, or effort offset ")
  print("5. VESSEL -> vessel ID ")
  print("Please ensure that latitude and depth in strata.limits match the following boundaries:")
  print("Latitude: 42-49 in 0.5 increments")
  print("Depth (meters): 55, 75, 100, 125, 155, 183, 200, 250, 300, 350, 400, 450, 500, 549, 600, 700, 800, 900, 1000, 1100, 1200, 1280")
  
  # Give information about necessary headers
  if(!all(c('BEST_DEPTH_M','BEST_LAT_DD',species,'YEAR','AREA_SWEPT_MSQ','VESSEL')%in%colnames(masterDat))){
    print("Warning: processData() terminated unsuccessfully.")
    print("Please ensure that masterDat has appropriate column names")
    stop()
  }
  
  # set up the generic data frame for this species
  Data = data.frame('PROJECT_CYCLE'=masterDat[,'YEAR'], 'BEST_DEPTH_M'=masterDat[,'BEST_DEPTH_M'], 'BEST_LAT_DD'=masterDat[,'BEST_LAT_DD'], 'HAUL_WT_KG'=masterDat[,which(dimnames(masterDat)[[2]]==species)], 'year'=as.factor(masterDat[,'YEAR']), 'effort'=masterDat[,'AREA_SWEPT_MSQ']*0.0001, 'VESSEL'=masterDat[,'VESSEL']) 
  if(Truncate>0){
    print(paste("Changing any observation with less than ",Truncate," kilograms to 0 kilograms",sep=""))
    Data[,'HAUL_WT_KG'] = ifelse( Data[,'HAUL_WT_KG']<Truncate, 0, Data[,'HAUL_WT_KG'] )
  }
  Data = cbind(Data, 'y'=Data[,'HAUL_WT_KG']) 
  Data = cbind(Data, 'strata'=apply(masterDat,1,strata.fn,Strata.df=strata.limits))
  Data = cbind(Data, 'isNonZeroTrawl'=ifelse(Data[,'y']>0,1,0))
  Data = cbind(Data, 'ones.vec'=rep(0,length(Data[,'y']))) # this is just for the 'ones-trick', inv Gaussian 
  Data = cbind(Data, 'logy3'=log(pi*2*(Data[,'y']^3))) # this is a constant for the invGaussian
  Data = cbind(Data, 'logy'=log((Data[,'y']))) # this is a constant for the LognormalECE2 model
  Data = cbind(Data, 'logeffort'=log(Data[,'effort']))
  Data = cbind(Data, 'effort2' = Data[,'effort']^2)
  Data = cbind(Data, 'logeffort2' = Data[,'logeffort']^2)  
  Data = cbind(Data, 'lfacty' = lfactorial(Data[,'y'])) 
  # Exclude tows from strata that are not included
  Exclude_NoStratum = which(is.na(Data$strata))
  print(paste("Excluded ",length(Exclude_NoStratum)," observations that were not assigned to any strata",sep=""))
  if(length(Exclude_NoStratum) > 0){
    print(paste("Observations that were not assigned to any strata are shown in 'Tows_outside_strata.csv'",sep=""))  
    write.table(Data[Exclude_NoStratum,],"Tows_outside_strata.csv",row.names=F,col.names=T,sep=",")
    Data = Data[-Exclude_NoStratum,]
  }
  
  # Exclude tows with some missing entry
  Exclude_Missing = which(apply(Data, MARGIN=1, FUN=function(Vec){any(is.na(Vec))}))
  print(paste("Excluded ",length(Exclude_Missing)," additional observations that had some missing data",sep=""))
  write.table(Data[Exclude_Missing,],"Tows_with_missing_data.csv",row.names=F,col.names=T,sep=",")
  if(length(Exclude_Missing) < 10 & length(Exclude_Missing) > 0) print(Data[Exclude_Missing,])
  if(length(Exclude_Missing) >= 10) print("Entries are not printed to the screen due to having 10 or more")
  if(length(Exclude_Missing) > 0) Data = Data[-Exclude_Missing,]
  
  # Count number of tows per strata and year
  TowsPerStrataYear = table(Data[,'strata'],Data[,'year'])
  write.table(TowsPerStrataYear,"Tows_Per_StrataYear.csv",row.names=F,col.names=T,sep=",")
  print(paste("Tows per strata and year are displayed below"))
  print(TowsPerStrataYear)
  
  # Count number of positive catches per strata and year
  EncountersPerStrataYear = table(Data[,'strata'],Data[,'year'],Data[,'isNonZeroTrawl'])[,,"1"]
  write.table(EncountersPerStrataYear,"Encounters_Per_StrataYear.csv",row.names=F,col.names=T,sep=",")
  print(paste("Encounters per strata and year are displayed below"))
  print(EncountersPerStrataYear)
  
  # Redefine variables that depend on year
  Data[,'year'] = factor(as.numeric(as.character(Data[,'year']))) # Record year factor in case years have been eliminated due to missing observations
  Data = data.frame(Data, 'vessel'=Letters[as.numeric(as.factor(as.character(Data[,'VESSEL'])))]) # Record year factor in case years have been eliminated due to missing observations
  # Define derived variables involving year
  Data = cbind(Data, 'strataYear'=factor(paste(Data[,'strata'],":",Data[,'year'],sep=""), levels=as.vector(outer(sort(unique(Data[,'strata'])),sort(unique(Data[,'year'])),FUN=paste,sep=":"))))
  Data = cbind(Data, 'vesselYear'=factor(paste(Data[,'vessel'],":",Data[,'year'],sep="")))
  # Attach and calculate other values that aren't in the data.frame
  #attach(Data)
  assign("Data", Data, envir = .GlobalEnv)
  #nonZeros = which(Data[,'isNonZeroTrawl']==TRUE)
  assign("nonZeros", which(Data[,'isNonZeroTrawl']==TRUE), envir = .GlobalEnv)
  # Number of elements
  #nS = length(unique(strata))
  assign("nS", length(unique(Data[,'strata'])), envir = .GlobalEnv)
  #nY = length(unique(year))
  assign("nY", length(unique(Data[,'year'])), envir = .GlobalEnv)
  #nNonZeros = length(nonZeros)
  assign("nNonZeros", length(nonZeros), envir = .GlobalEnv)
  #n = length(y)
  assign("n", length(Data[,'y']), envir = .GlobalEnv)
  #nSY = length(unique(strataYear)) 
  assign("nSY", nlevels(Data[,'strataYear']), envir = .GlobalEnv)
  #nVY = length(unique(vesselYear))
  assign("nVY", length(unique(Data[,'vesselYear'])), envir = .GlobalEnv)  
  #nV = length(unique(vesselYear))
  assign("nV", length(unique(Data[,'vessel'])), envir = .GlobalEnv)  
  # Diagonal matrix for the wishart / correlation model
  assign("R",diag(2), envir = .GlobalEnv)

  # If the covariates aren't in the R environment, create them
  if(!("X.bin" %in% ls(envir = .GlobalEnv))){
    assign("X.bin",matrix(NA, ncol=0, nrow=nrow(Data)), envir = .GlobalEnv)
  }else{
    if(length(Exclude_NoStratum) > 0) X.bin = X.bin[-Exclude_NoStratum]
    if(length(Exclude_Missing) > 0) X.bin = X.bin[-Exclude_Missing]
  }
  if(!("X.pos" %in% ls(envir = .GlobalEnv))){
    assign("X.pos",matrix(NA, ncol=0, nrow=nrow(Data)), envir = .GlobalEnv)
  }else{
    if(length(Exclude_NoStratum) > 0) X.pos = X.pos[-Exclude_NoStratum]
    if(length(Exclude_Missing) > 0) X.pos = X.pos[-Exclude_Missing]
  }
  nX.pos = ncol(X.pos)
  nX.binomial = ncol(X.bin)
  assign("nX.binomial",nX.binomial, envir = .GlobalEnv)
  assign("nX.pos",nX.pos, envir = .GlobalEnv)
          
  assign("log2pi",log(2*pi), envir = .GlobalEnv) # this is a constant for the LognormalECE2 model  
}