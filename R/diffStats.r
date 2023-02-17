#' Statistics on differentialy expressed regions.
#' 
#' @description 
#' `diffStats()` computes the percentage of rejected regions / bases.
#'
#' @param dds The `DESeqDataSet` object returned by [DiffSegR::dea()].
#' 
#' @export
diffStats <- function(dds){
  
  rejectedHypotheses <- GenomicRanges::mcols(dds)$rejectedHypotheses 
  
  if (is.null(rejectedHypotheses)) {
    stop("Please compute the differential expression analysis before.
         Results should be saved in `rejectedHypotheses`: TRUE -> 
         differentially expressed / FALSE -> not-differentially expressed.")
  }
  
  data.frame(
    features         = length(dds),
    featuresDiff     = length(dds[rejectedHypotheses,]),
    featuresDiffPerc = length(dds[rejectedHypotheses,])/length(dds),
    bases            = sum(GenomicRanges::width(dds)),
    basesDiff        = sum(GenomicRanges::width(dds[rejectedHypotheses,])),
    basesDiffPerc    = sum(GenomicRanges::width(dds[rejectedHypotheses,]))/
      sum(GenomicRanges::width(dds))
  )
}
