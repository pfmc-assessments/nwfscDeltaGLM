####################
# Function to calculate the correlations from the precision matrix
####################
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