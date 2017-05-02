#' Compute MLE index of abundance
#'
#' @param Data The data to be fit
#' @param Model The model object (as a list)
#' @param FileName The name of the file/directory
#' @param Folder Optional argument for named folder to store results (in current working directory)
#' @param Weights Defaults to "StrataAreas", other option is "Equal"
#' @param StrataTable Table of strata, required
#' @param Run Defaults to \code{TRUE}. Only run if there's no random effects and the distribution is either Gamma and Lognormal
#'
#' @return Returns a list of output, by year and by Strata:Year
#' @import R2jags
#' @import stats
#' @import utils
#' @export
#'
ComputeMleIndices = function(Data, Model, FileName, Folder=NA,
  Weights="StrataAreas", StrataTable, Run=TRUE){

  # Make folder
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")

  # Attach stuff -- listed by search()
  attach(Model$BUGSoutput$sims.list)
  #attach(Data)
  modelStructure = Model$modelStructure
  Dist = Model$likelihood

  # Estimate marginal means
  Mat = matrix(NA,nrow=length(levels(year)),ncol=length(levels(strata)))
  DataNew = data.frame(year=levels(year)[as.vector(row(Mat))], strata=levels(strata)[as.vector(col(Mat))])
  DataNew = data.frame(DataNew, logeffort=rep(log(1),nrow(DataNew)), vesselYear=rep(999,nrow(DataNew)), ones.vec=rep(log(1),nrow(DataNew)))   # I need to include ones.vec as zero for some reason

  # Only run if there's no random effects and the distribution is either Gamma and Lognormal
  if( (modelStructure$VesselYear.zeroTows%in%c("zero","fixed")) & (modelStructure$VesselYear.positiveTows%in%c("zero","fixed")) & (modelStructure$StrataYear.zeroTows%in%c("zero","fixed")) & (modelStructure$StrataYear.positiveTows%in%c("zero","fixed"))  & (Dist=="gamma" | Dist=="lognormal") & Run==TRUE){
    # Default formulae
    FormulaPres = " ~ 0 + factor(year)"
    if(nlevels(strata)>1) FormulaPres = paste(FormulaPres," + factor(strata)",sep="")
    FormulaPos = " ~ 0 + factor(year)"
    if(nlevels(strata)>1) FormulaPos = paste(FormulaPos," + factor(strata)",sep="")

    # Modified formulae
    if(modelStructure$StrataYear.zeroTows=="fixed" & nlevels(strata)>1){FormulaPres = paste(FormulaPres," + factor(strata):factor(year)",sep="")}
    #if(modelStructure$VesselYear.zeroTows=="fixed"){FormulaPres = paste(FormulaPres," + factor(vesselYear)",sep="")}
    #if(modelStructure$Catchability.zeroTows=="linear"){FormulaPres = paste(FormulaPres," + logeffort",sep="")}
    #  if(modelStructure$Catchability.zeroTows=="quadratic"){FormulaPres = paste(FormulaPres," + logeffort + logeffort2",sep="")}
    if(modelStructure$StrataYear.positiveTows=="fixed" & nlevels(strata)>1){FormulaPos = paste(FormulaPos," + factor(strata):factor(year)",sep="")}
    #if(modelStructure$VesselYear.positiveTows=="fixed"){FormulaPos = paste(FormulaPos," + factor(vesselYear)",sep="")}
    #if(modelStructure$Catchability.positiveTows=="linear"){FormulaPos = paste(FormulaPos," + logeffort",sep="")}
    #  if(modelStructure$Catchability.positiveTows=="quadratic"){FormulaPos = paste(FormulaPos," + logeffort + logeffort2",sep="")}

    #### Presence/absence
    GlmPres <- glm(as.formula(paste("ifelse(Data[,'HAUL_WT_KG']>0,1,0)",FormulaPres)), family=binomial, control=glm.control(epsilon=1e-8, maxit=1000, trace=FALSE))
    #### Positive catches
    #OffsetPos = list(ones.vec,logeffort)[[ifelse(modelStructure$Catchability.positiveTows=="one",2,1)]]
    if(Dist=="gamma") GlmPos <- glm(as.formula(paste("HAUL_WT_KG",FormulaPos)), family=Gamma(link="log"), offset=logeffort, subset=which(Data[,'HAUL_WT_KG']>0), control=glm.control(epsilon=1e-8, maxit=1000, trace=FALSE))
    if(Dist=="lognormal") GlmPos <- lm(as.formula(paste("log(HAUL_WT_KG)",FormulaPos)), offset=logeffort, subset=which(Data[,'HAUL_WT_KG']>0))

    # Calculate strata areas
    Area = matrix(NA,nrow=length(unique(StrataTable[,'year'])),ncol=length(unique(StrataTable[,'strata'])))
    for(YearI in 1:nrow(Area)){
      for(StratI in 1:ncol(Area)){
        AreaI = which(StrataTable[,'strataYear']==paste(levels(strata)[StratI],":",levels(year)[YearI],sep=""))
        if(Weights=="StrataAreas") Area[YearI,StratI] = StrataTable[AreaI,'Area_Hectares']
        if(Weights=="Equal") Area[YearI,StratI] = 1
      }}

    # Predict indices
    IndexPres = array(predict(GlmPres, newdata=DataNew, type="response"), dim=dim(Mat))
    IndexPos = array(predict(GlmPos, newdata=DataNew, type="response"), dim=dim(Mat))
    if(Dist=="lognormal") IndexPos = exp(IndexPos + summary(GlmPos)$sigma^2/2)  # Bias correction
    Index = Area * IndexPres * IndexPos

  }else{
    Area = IndexPres = IndexPos = Index = matrix(NA,nrow=length(unique(StrataTable[,'year'])),ncol=length(unique(StrataTable[,'strata'])))
  }

  # EJ's code
  #DglmData = data.frame(Catch=Data[,'HAUL_WT_KG']/Data[,'effort'], Year=Data[,'year'])
  #Dglm = dglm(DglmData, dist="gamma", J=TRUE)      # J=TRUE calculates jacknife estimates of CV

  # Return results
  Results1 = data.frame(year=levels(year)[as.vector(row(Mat))], strata=levels(strata)[as.vector(col(Mat))], Index=as.vector(Index), Pres=as.vector(IndexPres), Pos=as.vector(IndexPos))
  Results2 = data.frame(year=levels(year), Index=rowSums(Index), Pres=rowSums(IndexPres), Pos=rowSums(IndexPos))

  # Write and print output
  write.csv(Results1,file=paste(Folder,"/",FileName,"ResultsByYearAndStrata_MLE.csv",sep=""))
  write.csv(Results2,file=paste(Folder,"/",FileName,"ResultsByYear_MLE.csv",sep=""))

  # Detach stuff -- listed by search()
  #detach(Data)
  detach(Model$BUGSoutput$sims.list)

  return(list(byYearAndStrata = Results1, byYear = Results2))
}
