.onLoad = function(libname, pkgname){

options(stringsAsFactors=TRUE)
Letters = apply(MARGIN=1,FUN=paste,collapse="",expand.grid(letters,letters))

# assign to this environment to keep from overwriting user's workspace
assign("Letters", Letters, envir=.GlobalEnv)

}
