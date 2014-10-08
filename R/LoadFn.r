LoadFn <- function (file, ...) {
  ls.ext <- function(file) {
    local({
      base::load(file)
      base::ls()
    }) }
  
  base::load(file, .GlobalEnv, ...)
  ls.ext(file)
}