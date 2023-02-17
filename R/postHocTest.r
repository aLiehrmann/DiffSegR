#' postHocTest
#' 
#' Wrapper function around \code{sanssouci::posthocBySimes} 
#' (Blanchard et al., 2020, doi:10.1214/19-AOS1847) for multiple tests correction 
#' through post-hoc inference. 
#' 
#' @param SExp An object that inherits from the `SummarizedExperiment` class and
#' augmented with the results of the differential expression analysis. 
#' @param predicate A predicate function that returns a single TRUE or FALSE if
#' the hypothesis on which it is applied meets the conditions defined in the 
#' predicate. The criteria used in the predicate have to be defined for each
#' hypotheses in mcols(SExp).
#' @param alpha  A `Double`. The significance level of the (post-hoc) confidence 
#' bound on the type 1 error (see [sanssouci::posthocBySimes]).
#' @param tdpLowerBound A value for the minimum true positive rate 
#' on the returned set of rejected regions.
#' @param orderBy A `String`. The parameter on which regions are ordered before 
#' the post-hoc procedure. 
#' @param dichotomicSearch A `Boolean`. Use it to speed up the post-hoc procedure.
#' @param verbose A `Boolean`. Should all the operations performed be displayed ?
#' 
#' @return The `SummarizedExperiment` object returned by 
#' [DiffSegR::segmentation()] augmented with the metadata column 
#' `rejectedHypotheses`. The rejected regions belong to the largest set of 
#' selected regions with a true positive rate of at least `tdpLowerBound` 
#' on the `alpha` level.
#' 
#' @export
postHocTest <- function(
  SExp, 
  predicate, 
  alpha, 
  tdpLowerBound,
  orderBy          = "pvalue",
  verbose          = FALSE,
  dichotomicSearch = FALSE) {
  
  metadata            <- GenomicRanges::mcols(SExp)
  original_order      <- order(metadata[[orderBy]])
  metadata            <- metadata[original_order,]
  metadata            <- metadata[!is.na(metadata[[orderBy]]),]
  selected_hypotheses <- which(predicate(metadata))
  
  if (dichotomicSearch){
  #- assumption: tp_lower_bound is a decreasing function of p-val. ------------#
  #- therefore, we can use a dichotomic search. 
  	dicotomicSearch <- function(
  	  iLowerBound, 
  	  iUpperBound) {
  	  len <- iUpperBound - iLowerBound+1
  	  if (verbose) message(len)
  	  if (len>2){
  	    m <- round((iUpperBound+iLowerBound)/2)
  	    tp_lower_bound <- sanssouci::posthocBySimes(
  	      p      = metadata[[orderBy]], 
  	      select = selected_hypotheses[1:m], 
  	      alpha  = alpha,
  	      Rcpp = TRUE
  	    )
  	    tp_rate <- tp_lower_bound/m
  	    if (verbose) message(tp_rate)
  	    if (tp_rate == tdpLowerBound){ #- maybe accept threshold here -------#
  	      m
  	    } else if (tp_rate > tdpLowerBound) {
  	      dicotomicSearch(m+1,iUpperBound)
  	    } else {
  	      dicotomicSearch(iLowerBound,m-1)
  	    }
  	  } else {
        tp_lower_bound <- sanssouci::posthocBySimes(
  	      p      = metadata[[orderBy]], 
  	      select = selected_hypotheses[1:iLowerBound], 
  	      alpha  = alpha,
  	      Rcpp = TRUE
  	    )
        tp_upper_bound <- sanssouci::posthocBySimes(
  	      p      = metadata[[orderBy]], 
  	      select = selected_hypotheses[1:iUpperBound], 
  	      alpha  = alpha,
  	      Rcpp = TRUE
  	    )
        if (tp_upper_bound/iUpperBound>tdpLowerBound){
          iUpperBound
        } else if (tp_lower_bound/iLowerBound>tdpLowerBound){
          iLowerBound
        } else {
  	      iLowerBound-1
        }
  	  }
  	}
  	last_selected_hypothesis <- dicotomicSearch(1,length(selected_hypotheses))
  } else {
  	lower_bound_tp <- sapply(
			seq_along(selected_hypotheses), 
			function(x){
				sanssouci::posthocBySimes(
					p      = metadata[[orderBy]], 
					select = selected_hypotheses[1:x], 
					alpha  = alpha, 
					Rcpp   = TRUE
				)
			}
		)
		if (length(lower_bound_tp)>0){
			last_selected_hypothesis <- sum(
				lower_bound_tp/seq_along(selected_hypotheses)>tdpLowerBound
			)
		} else {
			last_selected_hypothesis <- c()
		}
  }
  tp <- sanssouci::posthocBySimes(
  	p      = metadata[[orderBy]], 
  	select = selected_hypotheses[1:last_selected_hypothesis], 
  	alpha  = alpha,
  	Rcpp = TRUE
  )
  rejected_hypotheses <- rep(FALSE, nrow(SExp))
  if (length(last_selected_hypothesis)>0){
  	rejected_hypotheses[selected_hypotheses[1:last_selected_hypothesis]] <- TRUE
  }
  GenomicRanges::mcols(SExp)$rejectedHypotheses <- rejected_hypotheses[order(original_order)]
  message(
    sum(GenomicRanges::mcols(SExp)$rejectedHypotheses), 
    " rejected hypotheses with at least ", 
    tp/(last_selected_hypothesis), 
    " of true positives."
  )
  SExp
}
