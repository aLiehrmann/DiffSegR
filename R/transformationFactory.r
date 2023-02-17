transformationFactory <- function(type) {
  fn <- switch(type,
    log2FoldChange = transformationFactory.log2FoldChange,
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