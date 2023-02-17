mergeFeatures <- function(
  SExp) {
  #- assumption: GRanges always sorted ---------------------------------------##
  features   <- as.data.frame(SummarizedExperiment::rowRanges(SExp))
  group      <- rep(NA, length(SExp))
  group[[1]] <- 0
  for (i in 2:nrow(features)){
    new_group     <- predicate(
      featureA = features[i-1,], 
      featureB = features[i,]
    )
    group[[i]] <- ifelse(new_group, group[[i-1]]+1, group[[i-1]])
    ##- cleaner option to try: -----------------------------------------------##
    ##- group[[i]] <- new_group; group <- cumsum(group) ----------------------##
  }
  SummarizedExperiment::mcols(SExp)$group <- group
  SExp
}

predicate <- function(featureA, featureB){
  if (featureA$strand != featureB$strand) {
    TRUE
  } else if (featureA$end+1 != featureB$start) {
    TRUE
  } else if (featureA$rejectedHypotheses != featureB$rejectedHypotheses) {
    TRUE
  } else if (
    featureB$rejectedHypotheses & 
    sign(featureA$log2FoldChangeMean) != 
    sign(featureB$log2FoldChangeMean)) {
    TRUE
  } else {
    FALSE
  }
}