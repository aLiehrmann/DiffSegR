#' denoise
#'
#' @description 
#' `denoise()` estimates the normalized denoised coverages and log2-FC 
#' per-base.
#' 
#' @param dds A `DESeqDataSet` object augmented with modelStart and modelEnd 
#' columns.
#' @param data The `List` object returned by [DiffSegR::loadData()].
#' 
#' @return A `DESeqDataSet` augmented with log2FCPerBase, log2RefPerBase, 
#' log2OtherPerBase columns.
#' 
#' @export
denoise <- function(
  dds,
  data){
  
  size_factors <- DESeq2::sizeFactors(dds)
  strands <- names(data$coverages)
  if (length(size_factors) != length(data$coverages[[strands[[1]]]])){
    stop("Please provide one library size factor by sample (sizeFactors).")
  }
  if (is.null(GenomicRanges::mcols(dds)$modelStart)){
    stop("No segment boundary found (modelStart, modelEnd).")
  }
  dico_strand <- list(
    "plus"  = "+",
    "minus" = "-",
    "all"   = "*"
  )
  referenceCondition <- data$referenceCondition
  i_reference        <- data$sampleInfo$condition == referenceCondition
  do.call(SummarizedExperiment::rbind, lapply(
    strands,
    function(strand){ 
      ##- list of RLE as matrix ----------------------------------------------##
      coverages <- do.call(cbind, lapply(
        data$coverages[[strand]],
        as.vector
      ))
      ##- compute normalized transformation & coverage -----------------------##
      scaled_coverages <- coverages%*%diag(1/size_factors)
      lfc_per_base <- transformationFactory(type = "log2FoldChange")(
        numerator   = scaled_coverages[,!i_reference, drop=FALSE],
        denominator = scaled_coverages[,i_reference, drop=FALSE]
      )
      log2_ref_per_base <- rowMeans(log2(
        scaled_coverages[,i_reference, drop=FALSE]+1
      ))
      log2_other_per_base <- rowMeans(log2(
        scaled_coverages[,!i_reference, drop=FALSE]+1
      ))
      ##- bounds of segments -------------------------------------------------##
      s_dds <- dds[GenomicRanges::strand(dds)==dico_strand[[strand]],]
      starts <- GenomicRanges::mcols(s_dds)$modelStart
      ends   <- GenomicRanges::mcols(s_dds)$modelEnd
      ##- estimate denoised transformation & coverage ------------------------##
      lfc_mean <- sapply(
        seq_along(starts),
        function(i) mean(lfc_per_base[starts[[i]]:ends[[i]]])
      )
      log2_other_mean <- sapply(
        seq_along(starts),
        function(i) mean(log2_other_per_base[starts[[i]]:ends[[i]]])
      )
      log2_ref_mean <- sapply(
        seq_along(starts),
        function(i) mean(log2_ref_per_base[starts[[i]]:ends[[i]]])
      )
      GenomicRanges::mcols(s_dds)$log2FoldChangeMean <- lfc_mean
      GenomicRanges::mcols(s_dds)$log2RefMean   <- log2_ref_mean
      GenomicRanges::mcols(s_dds)$log2OtherMean <- log2_other_mean
      s_dds
    }
  ))
}
