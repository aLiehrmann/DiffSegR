#' Summarize the differential transcription landscape
#' 
#' @description 
#' In brief, `segmentation()` (i) identifies segments with homogeneous log2-FC 
#' along the locus of interest and (ii) counts reads overlapping them.
#' 
#' In details, `segmentation()` (i) starts by calling [fpopw::Fpop_w()] on the 
#' log2-FC profile of both strands to retreive the segments boundaries 
#' (changepoints). The hyperparameter `alpha` specified by the user controls 
#' the number of returned segments. Corresponding genomic regions are stored as 
#' a GenomicRanges object. (ii) Finally, `segmentation()` calls 
#' [Rsubread::featureCounts()] which assigns to the regions the mapped reads 
#' from each replicate of each biological condition. By default a read is 
#' allowed to be assigned to more than one region if it is found to overlap 
#' with more than one region. [Rsubread::featureCounts()] ends by building a 
#' count matrix. `segmentation()` ends by returning the regions and associated
#' count matrix as `SummarizedExperiment` object.
#' 
#' @param data The `List` object returned by [DiffSegR::loadData()].
#' @param weightType A `String`. Select the type of weights associated to 
#' the log2-FC per-base:
#' \itemize{
#'   \item unweighted   : all observations have the same weight (=1) ; 
#'   \item zeroInflated : low counts have less weight.
#' }
#' @param modelSelectionType A `String`. Select the penalty used by FPOP:
#' \itemize{
#'   \item yao : Yao's penalty \code{alpha*sigma^2*log(n)}.
#' }
#' @param featureCountsType A `String`. Select how to summarize counts:
#' \itemize{
#'   \item fromBam      : from bam files using [Rsubread::featureCounts] ;
#'   \item fromCoverage : from coverage profiles.
#' }
#' @param compressed A `Boolean`. Indicate if the observations have to be
#' compressed (does not change the segmentation results and decreases
#' the running time).
#' @param outputDirectory A `String`. Path to the output directory.
#' @param alpha A `Double`. Segmentation hyperparameter used in Yao's penalty: 
#' \code{alpha*sigma^2*log(n)}. The number of changepoints returned by 
#' [fpopw::Fpop_w()] is a decreasing function of `alpha`.
#' @param segmentNeighborhood A `Boolean`. Indicate the weighted pDPA algorithm 
#' (Rigaill 2010 and 2015) has to be used to explore the space of segmentations.
#' @param Kmax An `Integer`. Segmentations with 1 to Kmax segments are recovered 
#' using segment neighborhood.
#' @param nbThreadsGridSearch An `Integer`. Number of threads used by the grid 
#' search procedure.
#' @param alphas A vector of `Double`. A series of alphas used by the grid 
#' search procedure.
#' @param gridSearch A `Boolean`. Indicate if a grid search has to be used 
#' to explore the space of segmentations.
#' @param strandSpecific An `Integer`. Passed on to [Rsubread::featureCounts()].
#' @param nbThreadsFeatureCounts An `Integer`. Passed on to 
#' [Rsubread::featureCounts()].
#' @param isPairedEnd Passed on to [Rsubread::featureCounts()].
#' @param read2pos An `Integer`. Passed on to [Rsubread::featureCounts()].
#' @param verbose A `Boolean`. Should all the operations performed be 
#' displayed ?
#' @param featureCountsOtherParams A `List`. Other paramters passed on to 
#' [Rsubread::featureCounts()].
#' @return A `SummarizedExperiment` object that contains region boundaries and 
#' associated count matrix.
#' 
#' @examples 
#' 
#' # Create a working directory for running the example.
#' working_directory <- "./DIFFSEGR_TEST"
#' dir.create(working_directory)
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
#' 
#' # Summarize the differential landscape.
#' SExp <- segmentation(
#'   data                   = data, 
#'   nbThreadsFeatureCounts = 1,
#'   outputDirectory        = working_directory
#' )
#' 
#' # In genomic order, first to fifth regions and ...
#' print(SummarizedExperiment::mcols(SExp)[1:5,])
#' # ... associated counts.
#' print(SummarizedExperiment::assay(SExp)[1:5,])
#' 
#' # delete working directory 
#' unlink(working_directory, recursive = TRUE)
#' 
#' @export
segmentation <- function(
  data,
  weightType         = "unweighted",
  modelSelectionType = "yao",
  featureCountsType  = "fromBam",
  
  ##- modelSelection parameters ----------------------------------------------##
  compressed          = TRUE,
  alpha               = 2,
  segmentNeighborhood = FALSE,
  Kmax                = NULL,
  verbose             = FALSE,
  nbThreadsGridSearch = 1,
  alphas              = NULL,
  gridSearch          = FALSE,
  outputDirectory     = '.',

  ##- featureCounts parameters -----------------------------------------------##
  nbThreadsFeatureCounts = 1,
  strandSpecific         = 1,
  read2pos               = NULL,
  isPairedEnd            = FALSE,
  featureCountsOtherParams = list()
  ) {
  
  if (segmentNeighborhood & is.null(Kmax)) {
    stop("Please provide an upper bound segmentation (Kmax) for the segment
      neighborhood procedure")
  }
  
  if (gridSearch & is.null(alphas)) {
    stop("Please provide a vector of alphas for the grid search procedure.")
  }
  
  referenceCondition <- data$referenceCondition
  ##- build your segmentation model ------------------------------------------##
  weightHeuristic <- weightFactory(type=weightType)
  modelSelection  <- modelSelectionFactory(type=modelSelectionType)
  
  ##- retrieve segment boundaries --------------------------------------------##
  start_time <- Sys.time() 
  model <- modelSelection(
    locus               = data$locus,
    coverages           = data$coverages,
    log2FoldChange      = data$log2FoldChange,
    sampleInfo          = data$sampleInfo,
    weightHeuristic     = weightHeuristic,
    compressed          = compressed,
    outputDirectory     = outputDirectory,
    alpha               = alpha,
    segmentNeighborhood = segmentNeighborhood,
    Kmax                = Kmax,
    verbose             = verbose,
    nbThreads           = nbThreadsGridSearch,
    alphas              = alphas,
    gridSearch          = gridSearch
  )
  
  ##- cast segments into regions ---------------------------------------------##
  features <- modelToGRanges(
    model = model,
    locus = data$locus
  )
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  message("finished in ", diff_time, " ", attr(diff_time, "units"))
  
  ##- Counting reads within region boundaries --------------------------------##
  featureCounts  <- featureCountsFactory(type=featureCountsType)
  start_time <- Sys.time()
  message("counting ...")
  feature_counts <- featureCounts(
    features       = features,
    coverages      = data$coverages,
    sampleInfo     = data$sampleInfo,
    nbThreads      = nbThreadsFeatureCounts,
    strandSpecific = strandSpecific,
    read2pos       = read2pos,
    isPairedEnd    = isPairedEnd,
    featureCountsOtherParams = featureCountsOtherParams
  )
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  message("finished in ", diff_time, " ", attr(diff_time, "units"))
  
  ##- return SummarizedExperiment --------------------------------------------##
  colData <- as.data.frame(lapply(data$sampleInfo, as.factor))  
  SExp <- SummarizedExperiment::SummarizedExperiment(
    assays    = feature_counts,
    rowRanges = features,
    colData   = colData
  )
  SummarizedExperiment::colData(SExp)$condition <- stats::relevel(
    x   = SummarizedExperiment::colData(SExp)$condition, 
    ref = referenceCondition
  )
  SExp
}

modelToGRanges <- function(
  model,
  locus) {
  d <- list(
    minus = "-",
    plus  = "+",
    all   = "*"
  )
  features <- do.call(rbind,lapply(
    names(model), 
    function(strand){
      model_starts <- c(1,model[[strand]]$changepointsVec[
        -length(model[[strand]]$changepointsVec)
      ]+1)
      model_ends   <- model[[strand]]$changepointsVec
      chrom_ends   <- model_ends+GenomicRanges::start(locus)-1
      chrom_starts <- model_starts+GenomicRanges::start(locus)-1
      model_df     <- data.frame(
        seqname     = locus@seqnames@values,
        start       = chrom_starts,
        end         = chrom_ends,
        strand      = d[[strand]],
        modelStart  = model_starts,
        modelEnd    = model_ends,
        modelMean   = model[[strand]]$means
      )
      model_df$featureId <- paste(
        model_df$seqname,
        model_df$start,
        model_df$end,
        strand,
        sep = "_"
      )
      model_df
    }
  ))
  GenomicRanges::makeGRangesFromDataFrame(
    features, 
    keep.extra.columns = TRUE
  )
}
