featureCountsFactory <- function(type) {
  fn <- switch(type,
    fromBam      = featureCountsFactory.fromBam,
    fromCoverage = featureCountsFactory.fromCoverage,
    ## user-defined couting method 
    eval(parse(text=type))
  )
  attr(fn, "type") <- type
  fn
}


featureCountsFactory.fromBam <- function(
  coverageDir,
  loci,
  sampleInfo,
  features,
  strandSpecific,
  isPairedEnd,
  nbThreadsByLocus = 1,
  nbThreads        = 1,
  featureCountsOtherParams = list()) {
  
  ## segments as Simple Annotation Format (SAF)
  annot <- data.frame(
    GeneID = features$featureId,
    Chr    = as.character(GenomeInfoDb::seqnames(features)),
    Start  = GenomicRanges::start(features),
    End    = GenomicRanges::end(features),
    Strand = as.character(GenomicRanges::strand(features))
  )
  
  featureCounts_args <- c(
    list(
      files             = sampleInfo$bam,
      allowMultiOverlap = TRUE,  
      annot.ext         = annot,
      strandSpecific    = strandSpecific,
      nthreads          = nbThreads,
      isPairedEnd       = isPairedEnd
    ),
    featureCountsOtherParams
  )

  ## counts by segments
  res <- do.call(
    Rsubread::featureCounts,
    featureCounts_args
  )

  colnames(res$counts) <- sampleInfo$sample 
  res$counts
}


featureCountsFactory.fromCoverage <- function(
  coverageDir,
  loci,
  sampleInfo,
  features,
  strandSpecific = NA,
  isPairedEnd = NA,
  nbThreadsByLocus = 1,
  nbThreads        = 1,
  featureCountsOtherParams = list()) {
  
  d <- list(
    `+` = "plus",
    `-` = "minus",
    `*` = "all"
  )

  do.call(rbind, customLapply(1:length(loci), function(i_locus) {
    
    current_locus <- loci[i_locus,]
    
    target_features <- features[features$parentLocus == current_locus$locusID,]  

    path_to_coverages <- file.path(
      coverageDir,
      paste0(
        current_locus$locusID, 
        "_", 
        sampleInfo$sample, 
        ".rds"
      )
    )
    
    coverage_by_sample <- loadRDS(pathToCoverages = path_to_coverages)

    coverage_by_strand <- formatCoverageList(
      sampleInfo       = sampleInfo,
      coverageBySample = coverage_by_sample
    )

    rm(coverage_by_sample)

    features_df <- as.data.frame(target_features)

    counts <- do.call(rbind, lapply(1:nrow(features_df), function(i_feature) {
      start   <- features_df$modelStart[[i_feature]]
      end     <- features_df$modelEnd[[i_feature]]
      strand  <- d[[features_df$strand[i_feature]]]
      samples <- names(coverage_by_strand[[strand]])
    
      ## sum reads overlapping i_feature in each sample
      sapply(
        samples,
        function(sample){
          as.integer(sum(coverage_by_strand[[strand]][[sample]][start:end]))
        }
      )
    })) 
  }, nbThreads = nbThreads %/% nbThreadsByLocus))
}