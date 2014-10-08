#test svn
updateDeltaGLMcode <- function (local = NULL, save = FALSE, revision = "newest", env="nwfscDeltaGLM", pos=2)
{
    if(!is.environment(env)) {  #if env is an environment simply wirte to that. If not, create a new environment
        if(!is.character(env)) stop("'env' must be an environment or a character vector. Use default if not sure to assign objects into global workspace\n")
        if(is.character(env)) {
            env <- attach(NULL,pos=pos,name=env)
        }
    }
    getwebnames <- function() {
        changes <- readLines("http://r-forge.r-project.org/scm/viewvc.php/pkg/nwfscDeltaGLM/R/?root=nwfscassmt")
        line <- changes[grep("Directory revision:", changes)+1]  #Most recent revision is on next line
        current_revision <- as.numeric(strsplit(strsplit(line, "revision=",fixed=T)[[1]][2],"\">")[[1]][1])  #"
        cat("current revision number:", current_revision, "\n")
        tmp <- readLines(paste("http://r-forge.r-project.org/scm/viewvc.php?view=rev&root=nwfscassmt&revision=",current_revision,sep=""))
        line <- tmp[grep("Date:",tmp)[1]+1]
        lastDate <- strsplit(strsplit(line,"<td>")[[1]][2]," ")[[1]]
        cat("most recent change:", lastDate[c(1:3,5:6)], "\n")
        if (revision == "newest") {
            webdir <- "http://r-forge.r-project.org/scm/viewvc.php/pkg/nwfscDeltaGLM/R/?root=nwfscassmt"
            revision <- current_revision
            revDate <- lastDate[c(1:3,5:6)]
        }else {
            stop("VERSIONS NOT WORKING")
            if (is.numeric(revision) && revision <= current_revision) {
                webdir <- paste("http://r4ss.googlecode.com/svn-history/r",
                  revision, "/trunk/", sep = "")
            }else{
                stop("'revision' input should either be 'newest', or an integer <",
                  current_revision)
            }
        }
        cat("getting file names from", webdir, "\n")
        lines <- readLines(webdir, warn = F)
        filenames <- lines[ union(grep("*\\.R\"", lines), grep("*\\.r\"", lines))]   #"
        filenames <- unlist(lapply(strsplit(filenames,"name=\""),function(x){x[2]}))
        filenames <- unlist(lapply(strsplit(filenames,"\""),function(x){x[1]}))
        return(list(filenames = filenames, revision = revision, current_revision=current_revision,revDate=revDate))
    }
    getwebfiles <- function(fileinfo) {
        filenames <- fileinfo$filenames
        revision <- fileinfo$revision
        current_revision <- fileinfo$current_revision
        revDate <- fileinfo$revDate
        if(revision!=current_revision) stop("Former revision unable to be read in\n")
        webdir <- "http://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/nwfscDeltaGLM/R/"   #GetAges.fn.R?revision=7&root=nwfscassmt
        n <- length(filenames)
        cat(n, "files found\n")
        if (save) {
            cat("saving all files to", local, "\n")
            cat("NWFSC Survey Code\n",file=paste(local,"RevisionInfo.txt",sep="/"))
            cat("This is revision",revision,"of",current_revision,"\n",file=paste(local,"RevisionInfo.txt",sep="/"),append=T)
            cat("From",revDate,"\n",file=paste(local,"RevisionInfo.txt",sep="/"),append=T)
        }
        for (i in 1:n) {
            webfile <- paste(webdir, filenames[i],"?revision=",revision,"&root=nwfscassmt", sep = "")
            if (filenames[i] == "updateSurveyCode.R")
                webfile <- paste(webdir, filenames[i],"?revision=",current_revision,"&root=nwfscassmt", sep = "")
            if (save) {
                localfile <- paste(local, filenames[i], sep = "/")
                temp <- readLines(webfile)
                writeLines(temp, localfile)
                cat("  writing ", localfile, "\n", sep = "")
            }
            else {
                cat("  sourcing ", filenames[i], "\n", sep = "")
                source(webfile,local=env)
            }
            flush.console()
        }
    }
    getlocalfiles <- function(local) {
        filenames <- dir(local, pattern = "*.R$")
        n <- length(filenames)
        cat(n, "files found in", local, "\n")
        for (i in 1:n) {
            cat("  sourcing ", filenames[i], "\n", sep = "")
            source(paste(local, filenames[i], sep = "/"),local=env)
            flush.console()
        }
    }
    if (is.null(local)) {
        fileinfo <- getwebnames()
        getwebfiles(fileinfo)
    }
    else {
        if (save) {
            fileinfo <- getwebnames()
            getwebfiles(fileinfo)
        }
        getlocalfiles(local)
    }
    cat("update complete.\n")
}
