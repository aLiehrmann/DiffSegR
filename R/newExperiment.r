#' Define a New Experiment
#' 
#' @description 
#' `newExperiment()` allows for loading various pieces of information related to RNA-Seq 
#' samples and the loci that will be analyzed.
#' 
#' @param sampleInfo A `String`. Path to the file with 
#' sample information. The file includes the following columns: 
#' \itemize{
#'   \item sample         : identifier for the sample ;
#'   \item condition      : identifier for the condition ;
#'   \item replicate      : identifier for the replicate ;
#'   \item bam            : path to the bam file ;
#'   \item coverage       : path to the coverage file (in rds format);
#'   \item isPairedEnd    : logical indicating if the library contain 
#'   paired-end reads or not (TRUE/FALSE) ;
#'   \item strandSpecific : integer indicating if a strand-specific analysis
#'   should be performed. ‘0’ (unstranded), ‘1’ (stranded) and ‘2’ (reversely
#'   stranded). WARNING: All samples must be consistently stranded or non-stranded.
#    Ensure uniformity across your dataset to avoid analysis errors.
#' } Alternatively, one may provide a `Data.frame` that contains the same columns.
#' @param referenceCondition A `String`. The reference condition of 
#' the RNA-Seq experiment. 
#' @param otherCondition A `String`. The alternative condition that will be
#' compare to the reference one.  
#' @param loci A `String`. Path to the file with loci information. The file 
#' includes the following columns: 
#' \itemize{
#'   \item seqid      : chromosome identifier (match one of the chromosome
#'   in the bam) for the target genomic region ; 
#'   \item chromStart : start position for the target genomic region ; 
#'   \item chromEnd   : end position for the target genomic region ; 
#'   \item locusID    : unique id for the target genomic region.
#' } Alternatively, one may provide a `Data.frame` that contains the same columns.
#' @param coverage A `String`. Path to the directory where the 
#' coverage files will be saved.
#' @param nbThreads An `Integer`. The overall number of threads used for the 
#' analysis.
#' @param nbThreadsByLocus An `Integer`. The number of threads dedicated to each
#' locus. nbThreadsByLocus cannot be greater than nbThreads. The number of loci 
#' analyzed in parallel is equal to the divisor of nbThreads by nbThreadsByLocus. 
#' @returns A `List`. The list containing loaded data.
#'
#' @export
newExperiment <- function(
  sampleInfo,
  referenceCondition,
  otherCondition,
  loci,
  coverage,
  nbThreads        = 1,
  nbThreadsByLocus = 1) {
  
  if (!dir.exists(coverage)) {
    stop("The provided directory path does not exist. Please provide a valid 
    directory path. This directory is essential for storing coverage files.")
  }

  if (nbThreads < nbThreadsByLocus) {
    stop("The total number of threads specified (nbThreads) is lower than the number 
    of threads allocated per locus (nbThreadsByLocus).")
  }

  if (is.character(sampleInfo)) {
    sampleInfo <- utils::read.table(
  	  file             = sampleInfo, 
  	  header           = TRUE, 
  	  stringsAsFactors = FALSE
    )
  }

  if (all(sampleInfo$condition != referenceCondition)) {
    stop("The ", referenceCondition, " condition provided does not match 
      any biological condition of the file with information on samples.")
  }

  if (all(sampleInfo$condition != otherCondition)) {
    stop(paste0("The ", otherCondition," condition does not match 
      any biological condition of the file with information on samples."))
  }

  ## 1:stranded, 2:reversely stranded, 0: unstranded
  stranded <- sampleInfo$strandSpecific %in% c(1,2) 
  if (!(all(stranded) | all(!stranded))) {
    stop("BAMs should be either all (reversely) stranded or all unstranded, 
    but not a mix of both.")
  } else {
    stranded <- stranded[[1]]
  }

  if (is.character(loci)) {
    loci <- utils::read.table(
  	  file             = loci, 
  	  header           = TRUE, 
  	  stringsAsFactors = FALSE
    )
    loci$chromStart <- as.integer(loci$chromStart)
    loci$chromEnd   <- as.integer(loci$chromEnd)
  }
  
  ## Cast interest locus as GRanges object.
  loci <- GenomicRanges::makeGRangesFromDataFrame(
    df                 = loci,
    keep.extra.columns = TRUE
  )

  ## Reducing loci allows for merging overlapping loci, but it drops the IDs 
  ## of each locus. 
  reduced_loci <- GenomicRanges::reduce(loci)
  
  ## Identify overlaps between the reduced loci and each locus it is composed 
  ## of, in order to merge their IDs as well.
  tab_overlaps <- GenomicRanges::findOverlaps(reduced_loci, loci)
  reduced_loci$locusID <- sapply(
    split(tab_overlaps, tab_overlaps@from), 
    function(x) {
      paste0(loci$locusID[x@to], collapse="_")
    }
  )

  list(
    sampleInfo         = sampleInfo,
    referenceCondition = referenceCondition,
    otherCondition     = otherCondition,
    loci               = loci,
    nbThreads          = nbThreads,
    nbThreadsByLocus   = nbThreadsByLocus,
    stranded           = stranded,
    coverageDir        = coverage
  )
}