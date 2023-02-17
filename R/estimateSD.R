#' Some robust SD estimators.
#' 
#' @description 
#' Robust SD estimators: MAD, weighted MAD, HALL.
#' 
#' @param y A vector of `Double`. Observations.
#' @param w A vector of `Double`. Weights associated to the observations.
#' @param method A `String`. The type of estimator: 
#' \itemize{
#'   \item "MAD"
#'   \item "wtd_MAD"
#'   \item "HALL"
#' }
#' @return SD estimation of y.
#' 
#' @export
SD <- function(
  y,
  w      = NULL, 
  method = 'MAD') {
  n <- length(y) 
  if (is.null(w)){
    w <- rep(1, n)
  }
  if(method == "MAD"){
    return(stats::mad(diff(y)/sqrt(2)))	
  }
  if(method == "wtd_MAD"){
    matrixStats::weightedMad(x=y, w=w)	
  }
  if(method == "HALL"){
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(y)
    mat[2, -n] = mat[2, -1]
    mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]   
    return(sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3)))
  }
}
