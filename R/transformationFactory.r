transformationFactory <- function(type) {
  fn <- switch(type,
    log2FoldChange     = transformationFactory.log2FoldChange,
    log2FoldChange.Rle = transformationFactory.log2FoldChange.Rle,
    avg.Rle            = transformationFactory.avg.Rle,
    ## user-defined transformation method 
    eval(parse(text=type))
  )
  attr(fn, "type") <- type
  fn
}

transformationFactory.log2FoldChange <- function(
  numerator, 
  denominator, 
  prior = 1) {
  rowMeans(log2(numerator+prior)) - rowMeans(log2(denominator+prior))
}


transformationFactory.log2FoldChange.Rle <- function(
  numerator, 
  denominator, 
  prior = 1) { 
  if (typeof(numerator) == "S4") {
    numerator <- list(numerator)
  }
  if (typeof(denominator) == "S4") {
    denominator <- list(denominator)
  }
  Reduce(f=`+`, x=lapply(numerator, function(cov) log2(cov+prior)))/length(numerator)-
  Reduce(f=`+`, x=lapply(denominator, function(cov) log2(cov+prior)))/length(denominator)
}


transformationFactory.avg.Rle <- function(
  coverage,
  prior = 1) {
  if (typeof(coverage) == "S4") {
    coverage <- list(coverage)
  }
  Reduce(f=`+`, x=lapply(coverage, function(x) log2(x+prior)))/length(coverage)
}