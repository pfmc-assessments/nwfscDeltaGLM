#' Load functions
#'
#' @param file file name to load
#' @param \\dots additional arguments
#' \describe{
#'   \item{directory}{The directory to save output in}
#'   \item{mods}{A list of one or more fitted delta-GLMM model objects}
#'   \item{McmcDiagnostics}{Boolean indicating whether or not to make the MCMCdiagnostics plots, defaults to FALSE.}
#' }
#' @export
#'
LoadFn <- function (file, ...) {
  ls.ext <- function(file) {
    local({
      base::load(file)
      base::ls()
    }) }

  base::load(file, .GlobalEnv, ...)
  ls.ext(file)
}
