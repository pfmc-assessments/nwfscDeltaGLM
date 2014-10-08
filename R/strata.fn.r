strata.fn <- function(x,Strata.df) {
  # function per A. Hicks, 5/5/2012
  tmpL <- as.numeric(x["BEST_LAT_DD"])>Strata.df$SLat & as.numeric(x["BEST_LAT_DD"])<=Strata.df$NLat
  tmpD <- as.numeric(x["BEST_DEPTH_M"])>Strata.df$MinDepth & as.numeric(x["BEST_DEPTH_M"])<=Strata.df$MaxDepth
  Char = as.character(Strata.df[tmpL&tmpD,"STRATA"]) 
  return(ifelse(length(Char)==0,NA,Char))
}