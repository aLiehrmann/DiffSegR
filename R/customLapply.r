customLapply <- function(X, FUN, nbThreads, info=NULL, ...) {
  
  if (nbThreads>1) {
    
    require(parallel)
    
    if(!is.null(info)) {
      message(info)
    }
    
    ## Shared memory (fork) is not available on windows.
    if (.Platform$OS.type == "unix") {
      cl <-parallel::makeForkCluster(nbThreads)
    } else {
      cl <-parallel::makePSOCKcluster(nbThreads)
    }
    
    on.exit(parallel::stopCluster(cl))
    parallel::parLapply(cl, X, FUN, ...)

  } else {
    lapply(X, FUN, ...)
  }
}