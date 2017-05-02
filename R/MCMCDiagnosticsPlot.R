#' Gelman-Rubin diagnostics. This function does MCMC diagnostics, producing several plots related to Gelman diagnostics, and optionally a plot for the Geweke statistics (though this doesn't converge sometimes)
#'
#' @param McmcList An MCMC list object of model parameters
#' @param parToMonitor Which parameter(s) of the fitted model object were monitored / estimated
#' @param FileName Species and model specific, where the outupt will be generated to
#' @param Folder Folder for output, defaults to working directory
#' @param Geweke Boolean, whether to also calculate Geweke diagnostics, defaults to FALSE
#'
#' @import grDevices
#' @export
#'
McmcDiagnosticsPlot = function(McmcList, parToMonitor, FileName="", Folder=NA, Geweke=FALSE){

  GelmanDiag = gelman.diag(McmcList[,parToMonitor])

  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")

  # Plot
  jpeg(
    filename = paste(Folder, FileName, "Gelman_points.jpg", sep = ""),
    width = 6,
    height = 6,
    units = "in",
    res = 200
  )
  par(mfrow=c(1,1))
  # all points below 1.05 shown as points, all above given names
  mycols = rep("black",length(GelmanDiag$psrf[,1]))
  mycols[which(GelmanDiag$psrf[,1] > 1.05)] = "white"
  plot(GelmanDiag$psrf[,1],xlab="psrf",ylab="Gelman: point estimate",main="",col=mycols,lwd=2)
  lines(c(-1000,1000),c(1.05,1.05),lwd=2,col="red")
  if(length(which(GelmanDiag$psrf[,1] > 1.05)) > 0) {
    text(seq(1,length(GelmanDiag$psrf[,1]))[which(GelmanDiag$psrf[,1] > 1.05)],y = GelmanDiag$psrf[which(GelmanDiag$psrf[,1] > 1.05),1],labels=names(GelmanDiag$psrf[,1])[which(GelmanDiag$psrf[,1] > 1.05)])
  }
  dev.off()

  # Plot histograms showing Gelman-Rubin scores
  jpeg(
    filename = paste(Folder, FileName, "Gelman_histograms.jpg", sep = ""),
    width = 5,
    height = 10,
    units = "in",
    res = 200
  )
  par(mfrow=c(2,1),mai=c(0.5,0.5,0.2,0.2))
  hist(GelmanDiag$psrf[,1],xlab="psrf",main="Gelman: point estimate")
  hist(GelmanDiag$psrf[,2],xlab="psrf",main="Gelman: upper C.I.")
  dev.off()

  # plot the geweke diagnostic   -- OFTEN DOESN"T CONVERGE
  if(Geweke==TRUE){
    GewekeDiag = geweke.diag(McmcList[,parToMonitor])
    jpeg(
      filename = paste(Folder, FileName, "Geweke.jpg", sep = ""),
      width = 8,
      height = 8,
      units = "in",
      res = 200
    )
    par(mfrow=c(3,1),mai=c(0.5,0.5,0.2,0.2))
    for(i in 1:3) {
      plot(unlist(GewekeDiag[[i]]), xlab=paste("Parameter, chain: ",i),ylab="Geweke Z-score",ylim=c(-3,3))
      lines(c(-1000,1000),c(1.96,1.96),col="red",lty=3)
      lines(c(-1000,1000),c(-1.96,-1.96),col="red",lty=3)
    }
    dev.off()
  }
}
