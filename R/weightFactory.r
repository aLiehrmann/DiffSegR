weightFactory <- function(type) {
  fn <- switch(type,
    unweighted   = weightFactory.unweighted,
    zeroInflated = weightFactory.zeroInflated,
    eval(parse(text=type))
  )
  attr(fn, "type") <- type
  fn
}

weightFactory.unweighted <- function(
  y,
  prior=1) {
  
  nb_rows <- nrow(y)
  rep(1, nb_rows)
}

weightFactory.zeroInflated <- function(
  y, 
  prior=1) {
  
  #- geometric mean of each row scaled by global geometric mean ------------##
  geomean_rows <- exp(rowMeans(log(y+prior)))-prior
  geomean_rows/(exp(mean(log(geomean_rows+prior)))-prior)
}
