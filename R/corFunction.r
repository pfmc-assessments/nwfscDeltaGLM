#' Function to calculate the correlations from the precision matrix
#'
#' @param McmcArray The array of MCMC output
#' @param Parameter The names of all parameters monitoried
#' @param this.names The name of the parameter of interest
#' @param Model The fitted model object
#'
#' @return Returns a list of output, by year and by Strata:Year
#' @export
#'
corFunction = function(McmcArray, Parameter,this.names, Model) {
  nchains = dim(McmcArray)[2]
  cors = matrix(NA, nrow = dim(McmcArray)[1], ncol = nchains)
  # Shortcut for calculating correlation of 2x2
  # Matrix A = [a b c d]
  # Ainv = 1/(ad -bc)*[d -b -c a], cor = -b/sqrt(d*a)
  for(i in 1:nchains) {
    thisMat = McmcArray[,i,grep(this.names,Model$Parameter)]
    cors[,i] = -thisMat[,3]/sqrt(thisMat[,1]*thisMat[,4])
  }
  return(cors)
}
