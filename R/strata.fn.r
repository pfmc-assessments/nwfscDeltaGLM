#' Assign latitude and longitude to strata
#'
#' @param x The raw data, with latitude in \code{BEST_LAT_DD} and longitude in \code{BEST_LON_DD} and depth in \code{BEST_DEPTH_M}
#' @param Strata.df Data frame of strata boundaries
#'
#' @return Char Character vector of strate assignments
#' @export
#'
strata.fn <- function(x, Strata.df) {
  # function per A. Hicks, 5/5/2012
  tmpL <- as.numeric(x["BEST_LAT_DD"]) > Strata.df$SLat & as.numeric(x["BEST_LAT_DD"]) <= Strata.df$NLat
  tmpD <- as.numeric(x["BEST_DEPTH_M"]) > Strata.df$MinDepth & as.numeric(x["BEST_DEPTH_M"]) <= Strata.df$MaxDepth
  Char <- as.character(Strata.df[tmpL & tmpD, "STRATA"])
  return(ifelse(length(Char) == 0, NA, Char))
}
