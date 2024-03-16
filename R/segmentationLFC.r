#' Segment the Per-base Log2-FC of All Loci
#' 
#' @description 
#' WARNING : [DiffSegR::coverage()] must be run before `segmentationLFC()`.
#'
#' In brief, `segmentationLFC()` identifies segments with homogeneous Log2-FC 
#' along the loci of interest. segmentationLFC() can operate on each locus in 
#' parallel.
#' 
#' In details, for each locus, `segmentationLFC()` calls [fpopw::Fpop_w()] on the 
#' Log2-FC profile (possibly stranded) to retreive the segments boundaries 
#' (changepoints). The reference condition specified by the user is set to 
#' the denominator of the Log2-FC. The hyperparameter `alpha` specified by 
#  the user controls the number of returned segments. Corresponding genomic 
#' regions are returned as `GRanges` object.
#' 
#' @param data The `List` object returned by [DiffSegR::newExperiment()].
#' @param modelSelectionType A `String`. Select the penalty used by FPOP:
#' \itemize{
#'   \item yao : Yao's penalty \code{alpha*sigma^2*log(n)}.
#' }
#' @param outputDirectory A `String`. Path to the output directory. On directory 
#' by locus will be create within `outputDirectory`.
#' @param alpha A `Double`. Segmentation hyperparameter used in Yao's penalty: 
#' \code{alpha*sigma^2*log(n)}. The number of changepoints returned by 
#' [fpopw::Fpop_w()] is a decreasing function of `alpha`.
#' @param segmentNeighborhood A `Logical`. Indicate the weighted pDPA algorithm 
#' (Rigaill 2010 and 2015) has to be used to explore the space of segmentations.
#' @param Kmax An `Integer`. Segmentations with 1 to Kmax segments are recovered 
#' using segment neighborhood.
#' @param alphas A vector of `Double`. A series of alphas used by the grid 
#' search procedure.
#' @param gridSearch A `Logical`. Indicate if a grid search has to be used 
#' to explore the space of segmentations.
#' @param verbose A `Logical`. Should all the operations performed be 
#' displayed ?
#' @return A `GRanges` object that contains the segment boundaries.
#'
#' @export
segmentationLFC <- function(
  data,
  modelSelectionType  = "yao",
  outputDirectory     = NULL,
  alpha               = 2,
  segmentNeighborhood = FALSE,
  Kmax                = NULL,
  alphas              = NULL,
  gridSearch          = FALSE,
  verbose             = TRUE) {

  if ((segmentNeighborhood | gridSearch)  & is.null(outputDirectory))
  {
    stop("Please provide an output directory for saving the models identified 
    during the segment neighborhood or grid search procedure.")
  }

  if (segmentNeighborhood & is.null(Kmax)) 
  {
    stop("Please provide an upper bound segmentation (Kmax) for the segment
      neighborhood procedure")
  }
  
  if (gridSearch & is.null(alphas)) 
  {
    stop("Please provide a vector of alphas for the grid search procedure.")
  }
  
  ## build your segmentation model
  modelSelection  <- modelSelectionFactory(type=modelSelectionType)
  
  target_samples <- data$sampleInfo[
    data$sampleInfo$condition %in% c(data$referenceCondition, data$otherCondition),]

  if (verbose & data$nbThreads %/% data$nbThreadsByLocus>1) {
    
    ## One message for all since the one inside the loop will not be displayed
    message("\n > Segmenting the log2-FC of all loci...")
  }

  start_time <- Sys.time() 
  
  #suppressWarnings(
  features <- do.call(c,customLapply(1:length(data$loci), function(i_locus) {
    
    current_locus <- data$loci[i_locus,] 

    if (verbose) {
      message(
        "\n > Segmenting the log2-FC of locus ",
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

    path_to_coverages <- sub(
      ".rds", 
      paste0("_", current_locus$locusID, ".rds"),
      target_samples$coverage
    )
    
    coverage_by_sample <- loadRDS(pathToCoverages = path_to_coverages)

    coverage_by_strand <- formatCoverageList(
      sampleInfo       = target_samples,
      coverageBySample = coverage_by_sample
    )

    rm(coverage_by_sample)

    lfc_by_strand <- buildLFC(
      coverageByStrand   = coverage_by_strand,
      referenceCondition = data$referenceCondition,
      sampleInfo         = target_samples
    )

    rm(coverage_by_strand)

    ## retrieve segment boundaries in the log2 fold change
    model <- modelSelection(
      locus               = current_locus,
      log2FoldChange      = lfc_by_strand,
      outputDirectory     = outputDirectory,
      alpha               = alpha,
      segmentNeighborhood = segmentNeighborhood,
      Kmax                = Kmax,
      alphas              = alphas,
      gridSearch          = gridSearch,
      verbose             = verbose
    )
    
    ## cast segments into genomic regions 
    features <- modelToGRanges(
      model = model,
      locus = current_locus
    )
    
    features
  }, nbThreads = data$nbThreads %/% data$nbThreadsByLocus))
  #)

  end_time  <- Sys.time()
  diff_time <- end_time - start_time
  if (verbose) message("\n > Finished in ", diff_time, " ", attr(diff_time, "units"))

  features
}

buildLFC <- function(
  coverageByStrand,  
  referenceCondition,
  sampleInfo) {

  i_reference <- referenceCondition == sampleInfo$condition
  
  lfc_by_strand  <- lapply(names(coverageByStrand), function(strand) {
      ## Do we really keep the RLE format here ?
      lfc <- transformationFactory(type = "log2FoldChange.Rle")(
        numerator   = coverageByStrand[[strand]][!i_reference],
        denominator = coverageByStrand[[strand]][i_reference],
        prior       = 1
      )
    }
  )

  names(lfc_by_strand) <- names(coverageByStrand)
  lfc_by_strand
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
    function(strand) {
      model_starts <- c(1,model[[strand]]$changepoints[-length(model[[strand]]$changepoints)]+1)
      model_ends   <- model[[strand]]$changepoints
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
      model_df$parentLocus <- locus$locusID
      model_df
    }
  ))
  GenomicRanges::makeGRangesFromDataFrame(
    features, 
    keep.extra.columns = TRUE
  )
}
