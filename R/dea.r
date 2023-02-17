#' Differential expression analysis (DEA)
#'
#' @description 
#' `dea` takes as input the recently built `SummarizedExperiment` object and calls 
#' `DESeq2` to test for each line, i.e. each region, the difference in expression
#' between the two compared biological conditions. The user can choose to control 
#' all the resulting p-values or only a subset of interest. In the first case 
#' a Benjamini–Hochberg procedure is used to control the False Discovery Rate 
#' (FDR). In the second case, the user specifies selection criteria for the 
#' regions, e.g. a threshold on the absolute log2-FC per-regions. Then, a 
#' post-hoc procedure is used to control the joint Family wise error Rate 
#' (jFWER) on the selected subset of regions.
#' 
#' @param data The `List` object returned by [DiffSegR::loadData()]. 
#' @param SExp The `SummarizedExperiment` object returned by 
#' [DiffSegR::segmentation()].
#' @param design A `formula` object or a `Matrix` object. Passed on to 
#' [DESeq2::DESeqDataSet].
#' @param sizeFactors A vector of `Double`. Sample-specific size factors.
#' @param significanceLevel A `Double`. The significance cutoff on q-values.
#' @param predicate A predicate function that returns a single TRUE or FALSE if
#' the region on which it is applied meets the conditions defined in the 
#' predicate. The criteria used in the predicate have to be defined for each
#' regions in mcols(SExp).
#' @param postHoc_significanceLevel  A `Double`. The significance level of the 
#' test procedure (see [sanssouci::posthocBySimes]).
#' @param postHoc_tdpLowerBound A `Double`. The minimum true positive 
#' proportion on the returned set of rejected regions.
#' @param orderBy A `String`. The parameter on which regions are sorted before 
#' the post-hoc procedure. 
#' @param verbose A `Boolean`. Should all the operations performed be displayed
#' (only post-hoc) ?
#' @param dichotomicSearch A `Boolean`. Use it to speed up the post-hoc procedure.
#' @return A`DESeqDataSet` object. 
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
#' # Build coverages and log2-FC per-base.
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
#' # Differential expression analysis. We control all the resulting 
#' # p-values using a Benjamini–Hochberg procedure at level 0.01.
#' dds <- dea(
#'   data              = data,
#'   SExp              = SExp,
#'   design            = ~ condition,
#'   sizeFactors       = rep(1,4),
#'   significanceLevel = 1e-2
#' )
#' 
#' # first to fifth regions 
#' print(SummarizedExperiment::mcols(dds)[1:5,])
#' 
#' # delete working directory 
#' unlink(working_directory, recursive = TRUE)
#' 
#' @export
dea <- function(
  data,
  SExp,
  design                    = ~ condition,
  sizeFactors               = NA,
  significanceLevel         = 0.05,
  predicate                 = NULL,
  postHoc_significanceLevel = 0.05,
  postHoc_tdpLowerBound     = 0.95,
  orderBy                   = "pvalue",
  verbose                   = FALSE,
  dichotomicSearch          = FALSE) {
 
  message("  > differential expression analysis ...")
  dds <- DESeq2::DESeqDataSet(
    se     = SExp, 
    design = design
  )
  if (any(is.na(sizeFactors))) {
    dds <- DESeq2::estimateSizeFactors(dds)
  } else {
    DESeq2::sizeFactors(dds) <- sizeFactors
  }
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)
  message("  >  multiple testing correction ...")
  GenomicRanges::mcols(dds) <- cbind(
    GenomicRanges::mcols(dds),
    DESeq2::results(dds)
  )
  message("  >  normalize denoised log2-FC and coverage ...")
  dds <- denoise(
    dds  = dds,
    data = data
  )
  if (is.null(predicate)){
    GenomicRanges::mcols(dds)$rejectedHypotheses <- ifelse(
      GenomicRanges::mcols(dds)$padj<significanceLevel,
      TRUE,
      FALSE
    )
    if (any(is.na(GenomicRanges::mcols(dds)$padj))) {
      GenomicRanges::mcols(dds[
        is.na(GenomicRanges::mcols(dds)$padj),
      ])$rejectedHypotheses <- FALSE 
    }
  } else {
    message("  >  postHoc selection ...")
    dds <- postHocTest(
      SExp             = dds, 
      predicate        = predicate,
      alpha            = postHoc_significanceLevel,
      tdpLowerBound    = postHoc_tdpLowerBound,
      orderBy          = orderBy,
      verbose          = verbose,
      dichotomicSearch = dichotomicSearch
    )
  }
  dds <- mergeFeatures(dds)
  dds
}
