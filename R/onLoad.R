.onLoad = function(libname, pkgname){
# Check installed packages, and install any missing
#ToInstall = setdiff( c("rjags","R2jags","runjags","superdiag","pscl","statmod","stats"), installed.packages() )
#if(length(ToInstall)>0) install.packages( ToInstall )

# load default data into workspace (THIS IS NOW DONE IN THE DESCRIPTION FILE
#require(rjags)
#require(R2jags)
#require(runjags)
#require(superdiag)
#require(pscl)
#require(statmod)
#require(stats)

options(stringsAsFactors=TRUE)
Letters = apply(MARGIN=1,FUN=paste,collapse="",expand.grid(letters,letters))

# assign to this environment to keep from overwriting user's workspace
assign("Letters", Letters, envir=.GlobalEnv)

}
