#' Calculate Coverage Profiles
#' 
#' @description 
#' `coverage()` calculates the coverage profiles for loci specified by the user, 
#' using data from BAM files.
#'
#' @param data The `List` object returned by [DiffSegR::newExperiment()].
#' @param subsettingBams A `Logical`. Determines whether the BAM files should 
#' be subset by the loci specified by the user. This can significantly improve 
#' computation time.
#' @param coverageType A `String`. Select how to compute the coverage profiles:
#' \itemize{
#'   \item fivePrime  : coverage profiles compute on 5' ends of reads ; 
#'   \item threePrime : coverage profiles compute on 3' ends of reads ;
#'   \item average    : coverageFactory coverage profile compute on 
#'   average of 5' & 3' ends of reads ;
#'   \item fullLength : coverage profiles compute on full length reads.
#' }
#' @param verbose A `Logical`. Should all the operations performed be displayed ?
#' @param featureCountsOtherParams A `List`. Other paramters passed on to 
#' [Rsubread::featureCounts()].
#'
#' @export 
coverage <- function(
  data,
  subsettingBams = FALSE,
  coverageType   = "average",
  verbose        = TRUE, 
  featureCountsOtherParams = list()) {
  
  ## choose coverage heuristic
  coverage_fn <- coverageFactory(type=coverageType)

  if (verbose & data$nbThreads %/% data$nbThreadsByLocus>1) {
    
    ## One message for all since the one inside the loop will not be displayed
    message("\n > Calculating coverage for all loci...")
  }

  start_time <- Sys.time() 
  
  tmp <- customLapply(1:length(data$loci), function(i_locus) {
    
    current_locus <- data$loci[i_locus,]  

    if (verbose) {
      message(
        "\n > Calculating coverage of locus ",
        as.character(GenomicRanges::seqnames(current_locus)),
        " ",
        GenomicRanges::start(current_locus),
        "-",
        GenomicRanges::end(current_locus),
        " (ID : ",
        current_locus$locusID,
        ") ..."
      )
    }

    coverage_by_sample <- coverage_fn(
      locus           = current_locus,
      stranded        = data$stranded,
      bams            = data$sampleInfo$bam,
      strandSpecific  = data$sampleInfo$strandSpecific,
      isPairedEnd     = data$sampleInfo$isPairedEnd,
      nbThreads       = data$nbThreadsByLocus,
      subsettingBams  = subsettingBams,
      tmpDirectory    = data$coverageDir,
      featureCountsOtherParams = featureCountsOtherParams
    )

    path_to_cov <- file.path(
      data$coverageDir,
      paste0(
        current_locus$locusID, 
        "_", 
        data$sampleInfo$sample, 
        ".rds"
      )
    )

    ## save coverages
    exportAsRDS(
      pathToCoverages  = path_to_cov,
      coverageBySample = coverage_by_sample
    )
  }, nbThreads = data$nbThreads %/% data$nbThreadsByLocus)

  end_time  <- Sys.time()
  diff_time <- end_time - start_time
  if (verbose) message("\n > Finished in ", diff_time, " ", attr(diff_time, "units"))
}