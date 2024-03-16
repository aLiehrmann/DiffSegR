modelSelectionFactory <- function(type) {
  fn <- switch(type,
    yao = modelSelectionFactory.yao,
    eval(parse(text=type))
    ## user-defined segmentation method 
  )
  attr(fn, "type") <- type
  fn
}

modelSelectionFactory.yao <- function(
  locus,
  log2FoldChange,
  outputDirectory     = NULL,
  alpha               = 2,
  segmentNeighborhood = FALSE,
  gridSearch          = FALSE,
  Kmax                = NA,
  alphas              = NULL,
  verbose             = TRUE) {
  
  selected_model <- list()
  for (strand in names(log2FoldChange)) {

    if (all(log2FoldChange[[strand]]@values==0)) {
      
      ## default model / no signal 
      selected_model[[strand]] <- list(
        K            = 1, 
        changepoints = length(log2FoldChange[[strand]]),
        penalty      = 0,
        means        = 0,
        loss         = 0
      )
    } else {    
      sigma <- matrixStats::weightedSd(
        log2FoldChange[[strand]]@values, 
        log2FoldChange[[strand]]@lengths
      )   

      ## exploring some models using segment neighborhood
      if (segmentNeighborhood) {
        all_models <- segmentNeighborhoodFn(
          y     = log2FoldChange[[strand]]@values,
          w     = log2FoldChange[[strand]]@lengths,
          sigma = sigma,
          Kmax  = Kmax
        )
      }

      ## exploring some models using grid search
      if (gridSearch) {
       	all_models <- gridSearchFn(
          y       = log2FoldChange[[strand]]@values,
          w       = log2FoldChange[[strand]]@lengths,
          alphas  = alphas,
          sigma   = sigma
        )
      } 

      ## estimate changepoints position for user-provided alpha
      selected_model[[strand]] <- wrapperFpopw(
        y       = log2FoldChange[[strand]]@values,
        w       = log2FoldChange[[strand]]@lengths,
        penalty = alpha*sigma^2*log(length(log2FoldChange[[strand]])) 
      )

      if (gridSearch | segmentNeighborhood) {
        
        ## add selected model to those find with gridSearch or segmentNeighborhood
        all_models <- rbind(
          all_models, 
          data.frame(
            K      = selected_model[[strand]]$K, 
            alphas = alpha
          )
        )

        ## plot and save the number of segments K as function of the 
        ## hyperparameter alpha
        KvsAlpha(
          locus           = locus,
          strand          = strand,
          K               = all_models$K,
          alphas          = all_models$alphas,
          targetAlpha     = alpha,
          targetK         = selected_model[[strand]]$K,
          outputDirectory = outputDirectory,
          verbose         = verbose
        )
      }

      ## uncompress changepoints
      uncompressed_pos <- cumsum(log2FoldChange[[strand]]@lengths)
      selected_model[[strand]]$changepoints <- uncompressed_pos[
        selected_model[[strand]]$changepoints
      ]
    }
  }

  selected_model
}

wrapperFpopw <- function(
  y, 
  w, 
  penalty) {

  ## estimate changepoints position
  fit <- fpopw::Fpop_w(
    x      = y, 
    w      = w, 
    lambda = penalty
  )
  
  ## segments start
  if (length(fit$t.est)==1) {
    model_starts <- 1
  } else {
    ## probably sufficient, remove test
    model_starts <- c(1, fit$t.est[-length(fit$t.est)]+1) 
  }

  ## segments end
  model_ends <- fit$t.est
  
  ## segments mean
  means      <- sapply(
    seq_along(model_starts),
    function(i){
      stats::weighted.mean(
        x = y[model_starts[[i]]:model_ends[[i]]], 
        w = w[model_starts[[i]]:model_ends[[i]]]
      )
    }
  )

  ## get mean at each position 
  y_means <- rep(means, times = diff(c(0, fit$t.est)))
  
  ## calculate loss 
  loss    <- sum(w*(y-y_means)^2)
  
  list(
    K               = fit$K,
    changepoints    = fit$t.est,
    penalty         = penalty,
    means           = means,
    loss            = loss
  )
}

segmentNeighborhoodFn <- function(
  y,
  w,
  sigma,
  Kmax) {

  J_est <- fpopw::Fpsn_w_nomemory(
    x    = y,
    w    = w,
    Kmax = Kmax+1
  )$J.est
  
  J_est    <- rev(J_est)
  J_star   <- J_est[[1]]
  K_star   <- length(J_est)
  pen_star <- 0

  ## introduce new candidate model at each iteration
  for (x in 2:length(J_est)) {

    ## find optimal model on [pen_star[[last]]; + inf] by comparing
    ## new model with pevious ones 
    last    <- length(J_star)
    new_pen <- (J_est[[x]]-J_star[[last]])/(K_star[[last]]-(length(J_est)-x+1))
    while(new_pen<=pen_star[[last]]) {
      last    <- last-1
      new_pen <- (J_est[[x]]-J_star[[last]])/(K_star[[last]]-(length(J_est)-x+1))
    }

    ## save last optimal model
    J_star    <- c(J_star[1:last],J_est[[x]])
    K_star    <- c(K_star[1:last],length(J_est)-x+1)
    pen_star <- c(
      pen_star[1:last],
      (J_star[[length(J_star)-1]]-J_star[[length(J_star)]])/(K_star[[length(K_star)]]-K_star[[length(K_star)-1]])
    )
  }

  pen_star <- pen_star[-1]
  K_star   <- K_star[-1]
  alphas   <- (pen_star + c(diff(pen_star)/2,1))/(sigma^2*log(sum(w)))
  data.frame(K=rev(K_star),alphas=rev(alphas))
}

gridSearchFn <- function(
  y,
  w,
  alphas,
  sigma) {
  ## estimate the number of segments for each user-provided alphas
  K  <- unlist(lapply(
    alphas, 
    function(alpha) {
      wrapperFpopw(
        y       = y,
        w       = w,
        penalty = alpha*sigma^2*log(sum(w))
      )$changepoints+1
    }
  ))
  data.frame(K,alphas)
}

KvsAlpha <- function(
  locus, 
  strand,
  K,
  alphas,
  targetAlpha,
  targetK,
  outputDirectory,
  verbose) {
  
  segments <- c() ## bait for devtools::check()
  g <- ggplot2::ggplot(
    data.frame(
      alphas   = alphas,
      segments = K
    ),
    ggplot2::aes(
      x = log2(alphas), 
      y = segments
    )
  ) + ggplot2::theme_bw()

  if (length(alphas)==1) {
    g <- g + ggplot2::geom_point()
  } else {
    g <- g + ggplot2::geom_line() + ggplot2::geom_point(size=1) 
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
  ggplot2::theme(
    text = ggplot2::element_text(size = 20)
  )

  current_output_dir <- file.path(outputDirectory, locus$locusID)
    if(!file.exists(current_output_dir)) {
      dir.create(current_output_dir)
    }

  grDevices::pdf(
    file.path(
      current_output_dir, 
      paste0(
        "all_models_",
        strand, 
        ".pdf"
      )
    ), 
    width  = 11, 
    height = 10
  )
  print(g)
  grDevices::dev.off()

  saveRDS(
    object = data.frame(K,alphas), 
    file   = file.path(
      current_output_dir, 
      paste0(
        "all_models_",
        strand, 
        ".rds"
      )
    )
  ) 
}
