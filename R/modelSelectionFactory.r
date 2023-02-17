modelSelectionFactory <- function(type) {
  fn <- switch(type,
    yao = modelSelectionFactory.yao,
    eval(parse(text=type))
  )
  attr(fn, "type") <- type
  fn
}

modelSelectionFactory.yao <- function(
  locus,
  coverages,
  log2FoldChange,
  sampleInfo,
  weightHeuristic,
  compressed          = TRUE,
  outputDirectory     = ".",
  alpha               = 2,
  segmentNeighborhood = FALSE,
  Kmax                = NA,
  verbose             = FALSE,
  gridSearch          = FALSE,
  nbThreads           = 1,
  alphas              = NA) {
  
  selected_model <- list()
  for (strand in names(coverages)){
    message(
    	"changepoint detection on ",
    	as.character(GenomeInfoDb::seqnames(locus)),
    	" ",
    	GenomicRanges::start(locus),
    	" ",
    	GenomicRanges::end(locus),
    	" ",
    	strand, 
    	" ..."
    )
    ##- weighting and compressing data ---------------------------------------##
    obs <- preprocessing(
      log2FoldChange  = log2FoldChange[[strand]],
      coverages       = coverages[[strand]],
      weightHeuristic = weightHeuristic,
      compressed      = compressed
    )
    sigma   <- matrixStats::weightedSd(obs$y, obs$w)
    n       <- sum(obs$compression_lengths)    
    all_models <- data.frame()
    ##- exploring models using segment neighborhood --------------------------##
    if (segmentNeighborhood) {
      message("  > exploring models using segment neighborhood ...")
      all_models <- segmentNeighborhoodFn(
        obs   = obs,
        sigma = sigma,
        n     = n,
        Kmax  = Kmax
      )
    }
    ##- exploring models using grid search -----------------------------------##
    if (gridSearch) {
      message("  > exploring models using grid search ...")
     	all_models <- gridSearchFn(
        obs       = obs,
        alphas    = alphas,
        nbThreads = 1,
        alpha     = alpha,
        sigma     = sigma,
        n         = n 
      )
    }  
    ##- model selection ------------------------------------------------------##
    message("  > model selection ...")
    penalty <- alpha*sigma^2*log(n)
    selected_model[[strand]] <- wrapperFpopw(
      y       = obs$y,
      w       = obs$w,
      penalty = penalty
    )
    all_models <- rbind(
      all_models, 
      data.frame(
        K      = selected_model[[strand]]$changepoints+1, 
        alphas = alpha
      )
    )
    ##- number of segments K as function of hyperparameter alpha  ------------##
    KvsAlpha(
      strand          = strand,
      K               = all_models$K,
      alphas          = all_models$alphas,
      targetAlpha     = alpha,
      targetK         = selected_model[[strand]]$changepoints+1,
      outputDirectory = outputDirectory
    )
    ##- uncompress changepoints ----------------------------------------------##
    if (compressed){
      uncompressed_pos <- cumsum(obs$compression_lengths)
      selected_model[[strand]]$changepointsVec <- uncompressed_pos[
        selected_model[[strand]]$changepointsVec
      ]
    }
    if (verbose) {
      message(
        "sigma: ",sigma,
        ', penalty: ',
        penalty,
        ', number of changepoints: ',
        selected_model$changepoints
      )
    }
  }
  selected_model
}

wrapperFpopw <- function(
    y, 
    w, 
    penalty) {
  fit <- fpopw::Fpop_w(
    x      = y, 
    w      = w, 
    lambda = penalty
  )
  if (length(fit$t.est)==1){
    model_starts <- 1
  } else {
    model_starts <- c(1, fit$t.est[-length(fit$t.est)]+1)
  }
  model_ends <- fit$t.est
  means      <- sapply(
    seq_along(model_starts),
    function(i){
      stats::weighted.mean(
        x = y[model_starts[[i]]:model_ends[[i]]], 
        w = w[model_starts[[i]]:model_ends[[i]]]
      )
    }
  )
  y_means <- rep(means, times = diff(c(0, fit$t.est)))
  loss    <- sum(w*(y-y_means)^2)
  list(
    changepoints    = fit$K-1,
    changepointsVec = fit$t.est,
    penalty         = penalty,
    means           = means,
    loss            = loss
  )
}

preprocessing <- function(
  log2FoldChange,
  coverages, 
  weightHeuristic, 
  priorCounts, 
  compressed) {

  coverages <- do.call(cbind,lapply(
    coverages, 
    as.vector
  ))
  w <- weightHeuristic(
    y     = coverages,
    prior = priorCounts
  )
  log2FoldChange <- as.vector(log2FoldChange)
  if (compressed) {
    breaks  <- which(diff(log2FoldChange)!=0 | diff(w)!=0)
    starts  <- c(1,breaks+1)
    ends    <- c(breaks, length(log2FoldChange))
    compression_lengths <- ends-starts+1
    list(
      y                   = log2FoldChange[starts],
      w                   = w[starts]*compression_lengths,
      compression_lengths = compression_lengths
    )
  } else {
    list(
      y                   = as.vector(log2FoldChange),
      w                   = w,
      compression_lengths = rep(1,length(log2FoldChange))
    )
  }
}

segmentNeighborhoodFn <- function(
  obs,
  sigma,
  n,
  Kmax) {
  
  res_sn <- fpopw::Fpsn_w_nomemory(
    x    = obs$y,
    w    = obs$w,
    Kmax = Kmax
  )
  J_est    <- res_sn$J.est
  pen_star <- sapply(
	  seq_along(J_est)[-length(J_est)], 
	  function(i){
	    max(sapply(
	      (i+1):length(J_est), 
	      function(j){
	        (J_est[j]-J_est[i])/(i-j)
	      }
	    ))
	  }
  )
  K <- unique(sapply(
    seq_along(pen_star), 
    function(i) which.min(pen_star[1:i])
  ))
  alphas <- pen_star[K]/(sigma^2*log(n))
  data.frame(K,alphas)
}

gridSearchFn <- function(
  strand,
  obs,
  alphas,
  alpha,
  sigma,
  n,
  nbThreads = 1) {
  cl <- parallel::makeCluster(nbThreads) #- Is it as fast as mclapply ? -------#
  K  <- unlist(parallel::parLapply(
    cl,
    alphas, 
    function(alpha){
      wrapperFpopw(
        y       = obs$y,
        w       = obs$w,
        penalty = alpha*sigma^2*log(n)
      )$changepoints+1
    }#, mc.cores=nbThreads
  ))
  parallel::stopCluster(cl)
  data.frame(K,alphas)
}

KvsAlpha <- function(
  strand,
  K,
  alphas,
  targetAlpha,
  targetK,
  outputDirectory = "."
  ){
  segments <- c() #- bait for devtools::check() -------------------------------#
  g <- ggplot2::ggplot(
    data.frame(
      alphas   = alphas,
      segments = K
    ),
    ggplot2::aes(
      x = log2(alphas), 
      y = segments
    )
  )+
  ggplot2::theme_bw()

  if (length(alphas)==1) {
    g <- g + ggplot2::geom_point()
  } else {
    g <- g + 
    ggplot2::geom_line() +
    ggplot2::geom_point(size=1) 
  }
  
  g <- g + ggplot2::geom_vline(
    xintercept = log2(targetAlpha), 
    color      = "red"
  )+
  ggplot2::geom_hline(
    yintercept = targetK, 
    color      = "blue"
  )+
  ggplot2::xlab("log2(alpha)")+
  ggplot2::ylab("number of segments (K)")+
  #ggplot2::scale_x_continuous(tr="log2")+
  ggplot2::theme(
    text = ggplot2::element_text(size = 20)
  ) 
  grDevices::pdf(
    file.path(
      outputDirectory, 
      paste0("all_models_", strand, ".pdf")
    ), 
    width  = 11, 
    height = 10
  )
  print(g)
  grDevices::dev.off()
  message(paste0(
    "  > report available in ",
    file.path(
      outputDirectory, 
      paste0("all_models_", strand, ".pdf")
    )
  ))
  saveRDS(
    data.frame(K,alphas), 
    file.path(
      outputDirectory,
      paste0("all_models_", strand,".rds")
  )) 
}
