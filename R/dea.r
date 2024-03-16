#' Differential Expression Analysis (DEA)
#'
#' @description 
#' `dea` takes as input the recently built `SummarizedExperiment` object and calls 
#' `DESeq2` to test for each line, i.e. each region, the difference in expression
#' between the two compared biological conditions. The user can choose to control 
#' all the resulting p-values or only a subset of interest. In the first case 
#' a Benjamini–Hochberg procedure is used to control the False Discovery Rate 
#' (FDR). In the second case, the user specifies selection criteria for the 
#' regions, e.g. a threshold on the absolute Log2-FC per-regions. Then, a 
#' post-hoc procedure is used to control the joint Family wise error Rate 
#' (jFWER) on the selected subset of regions.
#' 
#' @param SExp The `SummarizedExperiment` object returned by 
#' [DiffSegR::counting()].
#' @param design A `formula` object or a `Matrix` object. Passed on to 
#' [DESeq2::DESeqDataSet].
#' @param sizeFactors A vector of `Double`. Sample-specific size factors.
#' @param significanceLevel A `Double`. The significance cutoff on q-values.
#' @param predicate A predicate function that returns a single TRUE or FALSE if
#' the region on which it is applied meets the conditions defined in the 
#' predicate. The criteria used in the predicate have to be defined for each
#' regions in mcols(SExp).
#' @param postHoc_significanceLevel  A `Double`. The significance level of the 
#' test procedure (See [sanssouci::posthocBySimes]).
#' @param postHoc_tdpLowerBound A `Double`. The minimum true positive 
#' proportion on the returned set of rejected regions.
#' @param verbose A `Logical`. Should all the operations performed be displayed ?
#' @return A`DESeqDataSet` object augmented with the metadata column `DER`.
#' `DER` is set to `TRUE` if the region is called differentially expressed. 
#' 
#' @export
dea <- function(
  SExp,
  design                    = ~ condition,
  sizeFactors               = NA,
  significanceLevel         = 0.05,
  predicate                 = NULL,
  postHoc_significanceLevel = 0.05,
  postHoc_tdpLowerBound     = 0.95,
  verbose                   = TRUE) {

  ## differential expression analysis with DESEQ2 
  ## (See https://bioconductor.org/packages/release/bioc/html/DESeq2.html) 
  if (verbose) message("\n > differential expression analysis ...")
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
  if (verbose) message("\n > multiple testing correction ...")
  GenomicRanges::mcols(dds) <- cbind(
    GenomicRanges::mcols(dds),
    DESeq2::results(dds)
  )
  
  if (is.null(predicate)) {## user-defined predicate for posthoc ?
    
    ## this columns was previously name `rejectedHypotheses`, now `DER`
    GenomicRanges::mcols(dds)$DER <- ifelse(
      GenomicRanges::mcols(dds)$padj<significanceLevel,
      TRUE,
      FALSE
    )
    if (any(is.na(GenomicRanges::mcols(dds)$padj))) {
      
      ## in case of NA values, DER=FALSE
      GenomicRanges::mcols(dds[
        is.na(GenomicRanges::mcols(dds)$padj),
      ])$DER <- FALSE 
    }
  } else {
    
    if (verbose) message("\n > postHoc selection ...")
    dds <- postHocTest(
      SExp             = dds, 
      predicate        = predicate,
      alpha            = postHoc_significanceLevel,
      tdpLowerBound    = postHoc_tdpLowerBound,
      verbose          = verbose
    )
  }
  
  dds <- mergeFeatures(dds)
  dds
}


mergeFeatures <- function(SExp) {
  
  ## assumption: GRanges always sorted
  features   <- as.data.frame(SummarizedExperiment::rowRanges(SExp))
  group      <- rep(NA, length(SExp))
  group[[1]] <- 0
  
  for (i in 2:nrow(features)) {
    new_group     <- newGroupRule(
      featureA = features[i-1,], 
      featureB = features[i,]
    )
    group[[i]] <- ifelse(new_group, group[[i-1]]+1, group[[i-1]])
    ## cleaner option to try:
    ## group[[i]] <- new_group; group <- cumsum(group)
  }
  
  SummarizedExperiment::mcols(SExp)$group <- group
  SExp
}


newGroupRule <- function(featureA, featureB) {

  if (featureA$seqnames != featureB$seqnames) {
    TRUE
  } else if (featureA$strand != featureB$strand) {
    TRUE
  } else if (featureA$end+1 != featureB$start) {
    TRUE
  } else if (featureA$DER != featureB$DER) {
    TRUE
  } else if (
    featureB$DER & 
    sign(featureA$log2FoldChange) != 
    sign(featureB$log2FoldChange)) {
    TRUE
  } else {
    FALSE
  }
}
