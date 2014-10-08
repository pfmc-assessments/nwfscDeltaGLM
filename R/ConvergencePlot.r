##########
# Trace plots for each parameters
##########
ConvergencePlot = function(McmcArray, maxDims=8, parToMonitor, parnames, Nkeep=5000, FileName, Type="Trace", Folder=NA, Model){
  
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  Nparam = length(parToMonitor)
  Nplots = ceiling( Nparam / maxDims^2 )
  KeepSet = seq(1,dim(McmcArray)[1],length=min(dim(McmcArray)[1],Nkeep))
  
  if(Type%in%c("Trace","ACF")){
    for(PlotI in 1:Nplots){
      ParSet = (PlotI-1)*maxDims^2 + 1:min(Nparam-(PlotI-1)*maxDims^2,maxDims^2)
      Ncol=ceiling(sqrt(length(ParSet)))
      Nrow = ceiling(length(ParSet)/Ncol)
      jpeg(paste(Folder,FileName,"",Type,"_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,units="in",res=150)
      par(mfrow=c(Nrow,Ncol), mar=c(0,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      for(ParI in 1:length(ParSet)){
        if(Type=="Trace"){
          matplot(McmcArray[KeepSet,,parToMonitor[ParSet[ParI]]], type="l", lty="solid", col=rainbow(dim(McmcArray)[2],alpha=0.4), main=parnames[ParSet[ParI]], xaxt="n", xlab="",ylab="", lwd=2)
        }
        if(Type=="ACF"){
          Acf = apply(McmcArray[KeepSet,,parToMonitor[ParSet[ParI]],drop=FALSE],MARGIN=2,FUN=function(Vec){acf(Vec,plot=FALSE)$acf})
          if(!any(is.na(Acf))){ 
            matplot(Acf, type="h", lty="solid", col=rainbow(dim(McmcArray)[2],alpha=0.4), main=parnames[ParSet[ParI]], xaxt="n", xlab="",ylab="", ylim=c(0,1), lwd=2)
          }else{
            plot.new()
          }
        } 
      } # ParI loop
      dev.off()
    } # PlotI loop
  } 
  if(Type=="VarianceDensity"){
    ParSet = grep("sigma",parnames)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,FileName,"Variance_density.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    for(ParI in 1:length(ParSet)){
      plot(x=-999,y=-999,,main=parnames[ParSet[ParI]],xlab="",ylab="",xlim=range(McmcArray[,,parToMonitor[ParSet[ParI]]]),ylim=c(0,10/diff(range(McmcArray[,,parToMonitor[ParSet[ParI]]])))) 
      for(ChainI in 1:dim(McmcArray)[2]){
        lines(density(McmcArray[,ChainI,parToMonitor[ParSet[ParI]]]),col=rainbow(dim(McmcArray)[2],alpha=0.4)[ChainI])
      }
    }
    dev.off()
  }
  
  if(Type=="CorrelationDensity"){
    ParSet = grep("Tau",parnames) # divide by 4 because we're only plotting correlation
    Ncol=ceiling(sqrt(length(ParSet)/4))
    Nrow = ceiling((length(ParSet)/4)/Ncol)
    corMcmc = array(NA, dim = c(dim(McmcArray)[1],dim(McmcArray)[2], length(ParSet)/4))
    corNames = c("strataYearTau","vesselYearTau","strataTau","yearTau")
    kept = 0
    keptNames = ""
    for(i in 1:4) {
      if(length(grep(corNames[i],parnames)) > 0) {
        kept = kept + 1
        corMcmc[,,kept] = corFunction(McmcArray, Model$Parameter, corNames[i], Model)
        keptNames[kept] = corNames[i]
      }    	
    }
    plotNames = c("Strata-Year","Vessel-Year","strata","year")
    jpeg(paste(Folder,FileName,"Correlation_density.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    # create a new subdataset from just the correlations
    for(par in 1:kept) {
      dens = density(corMcmc[,1,par],from=-1,to=1) # this is for scaling ylim
      plot(corMcmc[1,1,1], xlim=c(-1,1), ylim=c(0, 1.1*max(dens$y/sum(dens$y))), col = "white", main = plotNames[par], ylab = "Density",xlab = "correlation")
      for(chainI in 1:dim(McmcArray)[2]) {
        # for each chain, plot a separate density plot
        dens = density(corMcmc[,chainI,par],from=-1,to=1)
        lines(dens$x, dens$y/sum(dens$y),col=rainbow(dim(McmcArray)[2],alpha=0.4)[chainI])    		
      }
    }
    dev.off()
  }
  
}