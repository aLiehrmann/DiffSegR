featureCountsFactory <- function(type) {
  fn <- switch(type,
    fromBam      = featureCountsFactory.fromBam,
    fromCoverage = featureCountsFactory.fromCoverage,
    eval(parse(text=type))
  )
  attr(fn, "type") <- type
  fn
}

featureCountsFactory.fromBam <- function(
  coverages,
  features,
  sampleInfo,
  nbThreads      = 1,
  strandSpecific = 1,
  read2pos       = NULL,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
  
  annot <- data.frame(
    GeneID = features$featureId,
    Chr    = as.character(GenomeInfoDb::seqnames(features)),
    Start  = GenomicRanges::start(features),
    End    = GenomicRanges::end(features),
    Strand = as.character(GenomicRanges::strand(features))
  ) 
  res <- do.call(
    Rsubread::featureCounts,
    c(
      list(
        files             = sampleInfo$bam,
        allowMultiOverlap = TRUE,  
        annot.ext         = annot,
        strandSpecific    = strandSpecific,
        nthreads          = nbThreads,
        read2pos          = read2pos,
        isPairedEnd       = isPairedEnd
      ),
      featureCountsOtherParams
    )
  )
  colnames(res$counts) <- sampleInfo$sample 
  res$counts
}

featureCountsFactory.fromCoverage <- function(
  coverages,
  features,
  sampleInfo,
  nbThreads      = 1,
  strandSpecific = 1,
  read2pos       = NULL,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
  
  d <- list(
    `+` = "plus",
    `-` = "minus",
    `*` = "all"
  )
  features_df <- as.data.frame(features)
  cl <- parallel::makeCluster(nbThreads)
  all_covs <- do.call(rbind,parallel::parLapply(
    cl,
    1:nrow(features_df), 
    function(i_feature){
      start   <- features_df$modelStart[[i_feature]]
      end     <- features_df$modelEnd[[i_feature]]
      strand  <- d[[features_df$strand[i_feature]]]
      samples <- names(coverages[[strand]])
      
      ##- sum of overlapping reads in each sample ----------------------------##
      sapply(
        samples,
        function(sample){
          as.integer(sum(coverages[[strand]][[sample]][start:end]))
        }
      )
    }#, mc.cores = nbThreads
  ))
  parallel::stopCluster(cl)
  all_covs
}
