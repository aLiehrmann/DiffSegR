#' Multiple Tests Correction Through Post-hoc Inference
#' 
#' @description 
#' Wrapper function around \code{sanssouci::posthocBySimes} 
#' (Blanchard et al., 2020, doi:10.1214/19-AOS1847) for multiple tests correction 
#' through post-hoc inference. 
#' 
#' @param SExp An object that inherits from the `SummarizedExperiment` class and
#' augmented with the results of the differential expression analysis. 
#' @param predicate A predicate function that returns a single TRUE or FALSE if
#' the regions on which it is applied meets the conditions defined in the 
#' predicate. The criteria used in the predicate have to be defined for each
#' regions in mcols(SExp).
#' @param alpha  A `Double`. The significance level of the (post-hoc) confidence 
#' bound on the type 1 error (See [sanssouci::posthocBySimes]).
#' @param tdpLowerBound A value for the minimum true positive rate 
#' on the returned set of rejected regions.
#' @param verbose A `Logical`. Should all the operations performed be displayed ?
#' 
#' @return The `SExp` input object augmented with the metadata column 
#' `DER`. `DER` is set to `TRUE` if the region is called differentially expressed.
#' 
#' @export
postHocTest <- function(
  SExp, 
  predicate, 
  alpha, 
  tdpLowerBound,
  verbose          = TRUE) {
	
	## Authors would like to thank Pierre Neuvial for his help with the use of 
	## the sansSouci package.

	## Simes family of thresholds (See 
	## https://doi.org/10.1093/bioinformatics/btac693)
	thr_simes <- alpha*(1:length(SExp))/length(SExp)

	
	## genomic regions meeting the user-defined thresholds
	selected_hypotheses <- which(predicate(GenomicRanges::mcols(SExp)))
  
	if (length(selected_hypotheses)==0) {
		
		if (verbose) message("\n > No genomic regions have met the thresholds defined by the user.")
		GenomicRanges::mcols(SExp)$DER <- FALSE

	} else {
		
		## upper bound for the number of false discoveries among most
		## significant regions
		fp <- sanssouci:::curveMaxFP(
			p.values = GenomicRanges::mcols(SExp)$pvalue[selected_hypotheses], 
			thr      = thr_simes
		)

		## size of the largest set of genomic regions meeting the user-defined 
		## thresholds as well as a lower TDP bound that also satisfies a threshold 
		## set by the user
		nb_rejected <- sum(1-fp/(1:length(fp))>=tdpLowerBound)
	
		if (nb_rejected>0) {
			
			## largest pvalue to reject amongst selected genomic regions
			bound <- sort(GenomicRanges::mcols(SExp)$pvalue[selected_hypotheses])[nb_rejected]

			GenomicRanges::mcols(SExp)$DER <- FALSE
			GenomicRanges::mcols(SExp)$DER[selected_hypotheses] <- TRUE
			GenomicRanges::mcols(SExp)$DER[GenomicRanges::mcols(SExp)$pvalue>bound] <- FALSE

			## if several genomic regions on the bound reject only one of them
			x <- GenomicRanges::mcols(SExp)$pvalue==bound & GenomicRanges::mcols(SExp)$DER
			GenomicRanges::mcols(SExp)$DER[x] <- c(
				TRUE, 
				rep(FALSE, length(GenomicRanges::mcols(SExp)$DER[x])-1)
			)

			if (verbose) {
				message(
					"\n > ",
  	  		sum(GenomicRanges::mcols(SExp)$DER), 
  	  		" DERs with at least ", 
  	  		1-fp[[nb_rejected]]/nb_rejected, 
  	  		" of true positives."
  			)
			}
		
		} else {
			
			GenomicRanges::mcols(SExp)$DER <- FALSE
			if (verbose)	message("\n > 0 rejected hypotheses.")
		
		}
	}
	
	return(SExp)
}
