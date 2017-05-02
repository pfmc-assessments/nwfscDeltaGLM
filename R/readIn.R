#' Scan data
#'
#' @param ncol number of columns
#' @param nlines Number of lines
#' @param \\dots Additional arguments passed to the function
#'
#' @return data frame of data
#' @export
#'
readIn <- function(ncol, nlines, ...) {
  x <- matrix(scan(,"", quiet=TRUE, nlines=nlines), ncol=ncol, byrow=TRUE)
  x <- as.data.frame(x, ..., stringsAsFactors = FALSE)
  names(x) <- x[1,]
  x <- x[-1,]
  rownames(x) <- 1:nrow(x)
  for( i in 1:ncol(x)){
    x[,i] <- if(all(is.na(suppressWarnings(as.numeric(x[,i]))))){ x[,i] }else{ as.numeric(x[,i]) }
  }
  return(x)
}
