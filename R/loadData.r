#' Compute the differential transcription profile (log2-FC per-base)
#' 
#' @description 
#' `loadData()` loads : 
#'  * the sample information related to an RNA-seq experiment ;
#'  * the coverage profile overlapping the locus specified by the user 
#'  (computed from bam files) ;
#'  * the log2-FC per-base (computed from coverages).
#' 
#' @param sampleInfo A `String`. Path to the file with 
#' sample information. The file includes the following columns: 
#' \itemize{
#'   \item sample    : identifier for the sample ;
#'   \item condition : identifier for the condition ;
#'   \item replicate : identifier for the replicate ;
#'   \item bam       : path to the bam file ;
#'   \item coverage  : path to the coverage file (in rds format).
#' }
#' @param referenceCondition A `String`. Reference condition of 
#' the RNA-Seq experiment. The reference condition specified by the 
#' user is set to the denominator of the log2-FC.
#' @param locus A `List`. The list contains the coordinates for the target 
#' genomic region: 
#' \itemize{
#'   \item seqid      : a `String`. Chromosome identifier
#'   for the target genomic region ; 
#'   \item chromStart : an `Integer`. Start position for the target 
#'   genomic region ; 
#'   \item chromEnd   : an `Integer`. End position for the target 
#'   genomic region.
#' }
#' @param fromBam A `Boolean`. Indicate if coverage profiles have to be
#' computed from bam files or loaded from pre-computed rds files.
#' @param stranded A `Boolean`. Indicate if the reads are stranded or not.
#' @param strandSpecific Passed on to [Rsubread::featureCounts()].
#' @param isPairedEnd Passed on to [Rsubread::featureCounts()].
#' @param coverageType A `String`. Select how to compute the coverage profiles:
#' \itemize{
#'   \item fivePrime  : coverage profiles compute on 5' ends of reads ; 
#'   \item threePrime : coverage profiles compute on 3' ends of reads ;
#'   \item average    : coverageFactory coverage profile compute on 
#'   average of 5' & 3' ends of reads ;
#'   \item center     : coverage profiles compute on center position of reads ;
#'   \item fullLength : coverage profiles compute on full length reads.
#' }
#' @param readLength An `Integer`. The average length of reads from bam files.
#' @param nbThreads Passed on to [Rsubread::featureCounts()].
#' @param verbose A `Boolean`. Should all the operations performed be displayed ?
#' @param featureCountsOtherParams A `List`. Other paramters passed on to 
#' [Rsubread::featureCounts()].
#' @returns A `List`. The list contains loaded data:
#' \itemize{
#'   \item coverage profiles  : a `List` of coverages as `Rle`;
#'   \item log2-FC            : an `Rle` of `Double` ;
#'   \item stranded           : see [DiffSegR::loadData()] ; 
#'   \item locus              : the target genomic region as `GRanges`; 
#'   \item sampleInfo         : see [DiffSegR::loadData()] ; 
#'   \item referenceCondition : see [DiffSegR::loadData()] ; 
#' }
#' 
#' @examples 
#' # Create a working directory for running the example.
#' working_directory <- "./DIFFSEGR_TEST"
#' dir.create(working_directory, showWarnings = FALSE)
#' 
#' # Save sample information in a text file.
#' sample_info <- data.frame(
#'   sample    = c("pnp1_1_1", "pnp1_1_2", "wt_1", "wt_2"),
#'   condition = rep(c("pnp1_1", "wt"), each = 2),
#'   replicate = rep(1:2,2),
#'   bam       = sapply(
#'     c("pnp1_1_1_ChrC_71950_78500.bam", 
#'       "pnp1_1_2_ChrC_71950_78500.bam",
#'       "wt_1_ChrC_71950_78500.bam",
#'       "wt_2_ChrC_71950_78500.bam"
#'      ),
#'      function(bam) system.file("extdata", bam, package = "DiffSegR")
#'   ),
#'   coverage  = file.path(
#'     working_directory,
#'     paste0(c("pnp1_1_1", "pnp1_1_2", "wt_1", "wt_2"), ".rds")
#'   )
#' )
#' print(sample_info)
#' write.table(
#'   sample_info, 
#'   file.path(working_directory, "sample_info.txt")
#' )
#' 
#' # Build coverages and log2-FC.
#' data <- loadData(
#'   sampleInfo         = file.path(working_directory,"sample_info.txt"),
#'   locus        = list(
#'     seqid      = "ChrC", 
#'     chromStart = 71950, 
#'     chromEnd   = 78500
#'   ),
#'   referenceCondition = "wt",
#'   stranded           = TRUE,
#'   fromBam            = TRUE,
#'   nbThreads          = 1
#' )
#' print(data)
#' 
#' # delete working directory 
#' unlink(working_directory, recursive = TRUE)
#' 
#' @export
loadData <- function(
  sampleInfo,
  referenceCondition,
  locus,
  stranded        = TRUE,
  strandSpecific  = 1,
  fromBam         = TRUE,
  coverageType    = "average",
  readLength      = 1,
  isPairedEnd     = FALSE,
  nbThreads       = 1,
  verbose         = TRUE,
  featureCountsOtherParams = list()) {
    
  ##- read information on samples --------------------------------------------##
  sample_info <- utils::read.table(
  	file             = sampleInfo, 
  	header           = TRUE, 
  	stringsAsFactors = FALSE
  )
  
  if (length(unique(sample_info$condition)) != 2) {
    stop("Please provide two biological conditions.")
  }
  
  if (all(sample_info$condition != referenceCondition)) {
    stop("The reference biological condition provided does not match 
      any biological condition of the file with information on samples.")
  }
  
  ##- cast interest locus as GRanges object ----------------------------------##
  locus <- GenomicRanges::makeGRangesFromDataFrame(locus)
  
  ##- choose coverage heuristic ----------------------------------------------##
  coverage_fn <- coverageFactory(type=coverageType)
  
  ##- create coverage profiles from BAM files or load them from previous run -##
  start_time <- Sys.time()
  
  if (fromBam) { 
    ##- parallel evaluation of bam files -------------------------------------##
    if (verbose) message("building coverage profiles ...")
    coverage_by_sample <- coverage_fn(
      locus           = locus,
      stranded        = stranded,
      bams            = sample_info$bam,
      nbThreads       = nbThreads,
      strandSpecific  = strandSpecific,
      isPairedEnd     = isPairedEnd,
      featureCountsOtherParams = featureCountsOtherParams
    )
    ##- save coverage profiles -----------------------------------------------##
    exportAsRDS(
      pathToCoverages  = sample_info$coverage,
      coverageBySample = coverage_by_sample
    )
  } else {
    ##- load pre-computed coverage profiles ----------------------------------##
    coverage_by_sample <- loadRDS(
      pathToCoverages = sample_info$coverage, 
      verbose         = verbose
    )
  }
  
  ##- format coverage profiles list ------------------------------------------##
  coverage_by_strand <- formatCoverageList(
    sampleInfo       = sample_info,
    coverageBySample = coverage_by_sample
  )
  
  ##- compute log2-FC --------------------------------------------------------##
  if (verbose) message("building log2-FC per-base ...")
  lfc_by_strand <- buildLFC(
    coverageByStrand   = coverage_by_strand,
    referenceCondition = referenceCondition,
    sampleInfo         = sample_info
  )
  
  end_time  <- Sys.time()
  diff_time <- end_time - start_time
  message("finished in ", diff_time, " ", attr(diff_time, "units"))  
  
  list(
    coverages    	     = coverage_by_strand,
    log2FoldChange     = lfc_by_strand,
    locus        	     = locus,
    stranded     	     = stranded,
    sampleInfo 	       = sample_info,
    referenceCondition = referenceCondition
  )
}

exportAsRDS <- function(
  pathToCoverages,
  coverageBySample) {
  for (i_coverage in seq_along(pathToCoverages)) {
    saveRDS(coverageBySample[[i_coverage]], pathToCoverages[[i_coverage]])
  }
}

loadRDS <- function(
  pathToCoverages, 
  verbose = TRUE) {
  coverage_by_sample <- lapply(
    seq_along(pathToCoverages),
    function(i_coverage_path){
      if (verbose) message(
        "loading coverage profile from ",
        pathToCoverages[[i_coverage_path]],
        " ",
        i_coverage_path,
        "/",
        length(pathToCoverages),
        " ..."
      )
      readRDS(pathToCoverages[[i_coverage_path]])
    }
  )
  coverage_by_sample
}

buildLFC <- function(
  coverageByStrand,  
  referenceCondition,
  sampleInfo) {
  i_reference <- referenceCondition == sampleInfo$condition
  lfc_by_strand  <- lapply( 
    names(coverageByStrand),
    function(strand) {
      coverages <- do.call(cbind,lapply(
        coverageByStrand[[strand]], 
        as.vector
      ))
      lfc <- S4Vectors::Rle(
        transformationFactory(type = "log2FoldChange")(
          numerator   = coverages[,!i_reference, drop=FALSE],
          denominator = coverages[,i_reference, drop=FALSE],
          prior       = 1
        )
      )
      lfc
    }
  )
  names(lfc_by_strand) <- names(coverageByStrand)
  lfc_by_strand
}

formatCoverageList <- function(
  sampleInfo,
  coverageBySample) {
  
  samples                   <- sampleInfo$sample
  names(coverageBySample)   <- samples
  strands                   <- names(coverageBySample[[1]])
  coverage_by_strand        <- vector("list", length(strands))
  names(coverage_by_strand) <- strands 
  
  for (strand in strands){
    coverage_by_strand[[strand]]        <- vector("list", length(samples))
    names(coverage_by_strand[[strand]]) <- samples
    for (sample in samples){
      coverage_by_strand[[strand]][[sample]] <- 
        coverageBySample[[sample]][[strand]]
    }
  }
  coverage_by_strand
}
