PlotOffset = function(Data, BugsList, maxDims=8, FileName, Folder=NA){
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff
  attach(BugsList)
  #attach(Data)
  nonZeros = which(isNonZeroTrawl==TRUE)
  
  # Positive offset
  LogEffortRange = seq(min(logeffort),max(logeffort),length=1000)
  Nstrat = length(unique(strataYear[nonZeros]))
  Nplots = ceiling( Nstrat / maxDims^2 )
  for(PlotI in 1:Nplots){
    ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,"/",FileName,"Offset_Positive_by_strata_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=150,units="in")
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    for(StrataYearI in ParSet){
      StrataI = strata[which(strataYear==unique(strataYear)[StrataYearI])[1]]
      YearI = year[which(strataYear==unique(strataYear)[StrataYearI])[1]]
      Y = exp( median(cMx(Sdev)[,StrataI]) + median(Ydev[,YearI]) + median(SYdev[,StrataYearI]) + median(B.pos[,1])*LogEffortRange + median(B.pos[,2])*LogEffortRange )
      plot(x=exp(LogEffortRange),y=Y,ylab="",xlab="",main=unique(strataYear[nonZeros])[StrataI],ylim=c(0,max(Y)),xlim=c(0,max(exp(LogEffortRange))), type="l")
    }
    dev.off()
  }
  
  # Presence/Absence offset
  LogEffortRange = seq(min(logeffort),max(logeffort),length=1000)
  Nstrat = length(unique(strataYear[nonZeros]))
  Nplots = ceiling( Nstrat / maxDims^2 )
  for(PlotI in 1:Nplots){
    ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,"/",FileName,"Offset_Presence-Absence_by_strata_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=150,units="in")
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    for(StrataYearI in ParSet){
      StrataI = strata[which(strataYear==unique(strataYear)[StrataYearI])[1]]
      YearI = year[which(strataYear==unique(strataYear)[StrataYearI])[1]]
      Y = plogis( median(cMx(pSdev)[,StrataI]) + median(pYdev[,YearI]) + median(pSYdev[,StrataYearI]) + median(B.zero[,1])*LogEffortRange + median(B.zero[,2])*LogEffortRange )
      plot(x=exp(LogEffortRange),y=Y,ylab="",xlab="",main=unique(strataYear[nonZeros])[StrataI],ylim=c(0,1),xlim=c(0,max(exp(LogEffortRange))), type="l")
    }
    dev.off()
  }
  
  # Detach stuff
  #detach(Data)
  detach(BugsList)
}
