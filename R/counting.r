#' Quantifying Region Expression 
#'
#' @description 
#' In brief, `counting()` quantify expression of regions found by segmenting
#' the per-base Log2-FC of each locus. 
#'
#' In details `counting()` calls  [Rsubread::featureCounts()] which assigns to 
#' the regions the mapped reads from each replicate of both the reference and the
#' alternative biological condition. By default a read is allowed to be assigned 
#' to more than one region if it is found to overlap with more than one region. 
#' [Rsubread::featureCounts()] ends by building a count matrix. `counting()` then
#' returns regions and associated count matrix as `SummarizedExperiment` object.
#' Alternatively, users may opt to quantify the expression of each region 
#' directly from the coverage profiles.
#'
#' @param data The `List` object returned by [DiffSegR::newExperiment()].
#' @param features A `GRanges` object that contains the segment boundaries
#' (returned by [DiffSegR::segmentationLFC()]).
#' @param featureCountsType A `String`. Select how to summarize counts:
#' \itemize{
#'   \item fromBam      : from bam files using [Rsubread::featureCounts] ;
#'   \item fromCoverage : from coverage profiles.
#' }
#' @param verbose A `Boolean`. Should all the operations performed be 
#' displayed ?
#' @param featureCountsOtherParams A `List`. Other paramters passed on to 
#' [Rsubread::featureCounts()].
#' @export
counting <- function(
  data,
  features,
  featureCountsType        = "fromBam",
  featureCountsOtherParams = list(),
  verbose                  = TRUE) {
  ## counting reads within region boundaries
  featureCounts  <- featureCountsFactory(type=featureCountsType)

  target_samples <- data$sampleInfo[
    data$sampleInfo$condition %in% c(data$referenceCondition, data$otherCondition),]

  if (verbose) message("\n > Counting ...")
  
  start_time     <- Sys.time()

  feature_counts <- featureCounts(
    coverageDir              = data$coverageDir,
    loci                     = data$loci,
    sampleInfo               = target_samples,
    features                 = features,
    nbThreads                = data$nbThreads,
    nbThreadsByLocus         = data$nbThreadsByLocus,
    isPairedEnd              = target_samples$isPairedEnd,
    strandSpecific           = target_samples$strandSpecific,
    featureCountsOtherParams = featureCountsOtherParams
  )

  end_time  <- Sys.time()
  diff_time <- end_time - start_time
  if (verbose) message("\n > Finished in ", diff_time, " ", attr(diff_time, "units"))
  
  ## return SummarizedExperiment object
  colData <- as.data.frame(lapply(target_samples, as.factor))  
  feature_counts <- feature_counts[features$featureId,] #force the same order
  SExp    <- SummarizedExperiment::SummarizedExperiment(
    assays    = feature_counts,
    rowRanges = features,
    colData   = colData
  )

  SummarizedExperiment::colData(SExp)$condition <- stats::relevel(
    x   = SummarizedExperiment::colData(SExp)$condition, 
    ref = data$referenceCondition
  )
  
  SExp
}