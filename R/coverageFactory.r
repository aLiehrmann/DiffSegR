coverageFactory <- function(type) {
  fn <- switch(type,
    fivePrime   = coverageFactory.fivePrime,
    threePrime  = coverageFactory.threePrime,
    average     = coverageFactory.average,
    fullLength  = coverageFactory.fullLength,
    ## user-defined coverage calculation method 
    eval(parse(text=type)) 
  )
  attr(fn, "type") <- type
  fn
}

coverageFactory.fivePrime <- function(
  locus,
  bams,
  strandSpecific,
  isPairedEnd,
  verbose = TRUE) {
  
  param <- Rsamtools::ScanBamParam(which = locus)
  
  all_cov <- lapply(seq_along(bams), function(i) {
    if (verbose) cat(bams[[i]], "...\n")
 
    if (isPairedEnd[[i]]) {
      
      suppressWarnings(mates <- GenomicAlignments::readGAlignmentPairs(
        bams[[i]], 
        param      = param, 
        strandMode = strandSpecific[[i]]
      )) 
      gc()
      pairedEndCov(locus, mates, "5", strandSpecific[[i]])
            
    } else {
      
      reads <- GenomicAlignments::readGAlignments(
        bams[[i]], 
        param = param
      )
      gc()
      singleReadCov(locus, reads, "5", strandSpecific[[i]])
    }
  })
     
  names(all_cov) <- bams
  all_cov
}


coverageFactory.threePrime <- function(
  locus,
  bams,
  strandSpecific,
  isPairedEnd,
  verbose = TRUE) {
  
  param <- Rsamtools::ScanBamParam(which = locus)
  
  all_cov <- lapply(seq_along(bams), function(i) {
    if (verbose) cat(bams[[i]], "...\n")
 
    if (isPairedEnd[[i]]) {
      
      suppressWarnings(mates <- GenomicAlignments::readGAlignmentPairs(
        bams[[i]], 
        param      = param, 
        strandMode = strandSpecific[[i]]
      )) 
      gc()
      pairedEndCov(locus, mates, "3", strandSpecific[[i]])
            
    } else {
      
      reads <- GenomicAlignments::readGAlignments(
        bams[[i]], 
        param = param
      )
      gc()
      singleReadCov(locus, reads, "3", strandSpecific[[i]])
    }
  })
    
  names(all_cov) <- bams
  all_cov
} 


coverageFactory.average <- function(
  locus,
  bams,
  strandSpecific,
  isPairedEnd,
  verbose = TRUE) {
  
  param <- Rsamtools::ScanBamParam(which = locus)
  
  all_cov <- lapply(seq_along(bams), function(i) {
    if (verbose) cat(bams[[i]], "...\n")
    
    if (isPairedEnd[[i]]) {
      
      suppressWarnings(mates <- GenomicAlignments::readGAlignmentPairs(
        bams[[i]], 
        param      = param, 
        strandMode = strandSpecific[[i]]
      )) 
      gc()
      cov_ext <- list(
        threePrime = pairedEndCov(locus, mates, "3", strandSpecific[[i]]),
        fivePrime  = pairedEndCov(locus, mates, "5", strandSpecific[[i]])
      )
            
    } else {
      
      reads <- GenomicAlignments::readGAlignments(
        bams[[i]], 
        param = param
      )
      gc()
      cov_ext <- list(
        threePrime = singleReadCov(locus, reads, "3", strandSpecific[[i]]),
        fivePrime  = singleReadCov(locus, reads, "5", strandSpecific[[i]])
      )
    }
    
    cov_avg <- list()
    
    for (strand in names(cov_ext[[1]])) {
        cov_avg[[strand]] <- exp((log(cov_ext[["threePrime"]][[strand]]+1)+
          log(cov_ext[["fivePrime"]][[strand]]+1))/2)-1
    }
    
    cov_avg
  })
  
  names(all_cov) <- bams
  all_cov
} 

coverageFactory.fullLength <- function(
  locus,
  bams,
  strandSpecific,
  isPairedEnd,
  verbose = TRUE) {
  
  param <- Rsamtools::ScanBamParam(which = locus)
  
  all_cov <- lapply(seq_along(bams), function(i) {
    if (verbose) cat(bams[[i]], "...\n")
    
    if (isPairedEnd[[i]]) {
      
      suppressWarnings(reads <- GenomicAlignments::readGAlignmentPairs(
        bams[[i]], 
        param      = param, 
        strandMode = strandSpecific[[i]]
      )) 
      gc()
            
    } else {
      
      reads <- GenomicAlignments::readGAlignments(
        bams[[i]], 
        param = param
      )
      gc()
      if (strandSpecific[[i]] == 2) {
        reads <- BiocGenerics::invertStrand(reads)
      }
      
    }
    
    if (strandSpecific[[i]] == 1 | strandSpecific[[i]] == 2) {
      cov_by_strand <- lapply(c("+", "-"), function(s) {
        cov <- GenomicAlignments::coverage(reads[BiocGenerics::strand(reads)==s,])
        cov[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)]  
      }) 
      names(cov_by_strand) <- c("plus", "minus")
      cov_by_strand
    } else {
      list(all = GenomicAlignments::coverage(mates))
    }
  })
  names(all_cov) <- bams
  all_cov
} 



pairedEndCov <- function(locus, mates, typeExt, strandSpecific) {
  # we retrieve the first and second reads as well as the strand expression information
  first    <- GenomicAlignments::first(mates, real.strand=TRUE)
  last     <- GenomicAlignments::last(mates,  real.strand=TRUE)
      
  # we reduce the reads to their 5' or 3' ends while taking the strand into account
  first_ext <- GenomicRanges::resize(GenomicRanges::GRanges(first), width=1, fix=ifelse(typeExt=="5", "start", "end"))
  last_ext  <- GenomicRanges::resize(GenomicRanges::GRanges(last),  width=1, fix=ifelse(typeExt=="5", "start", "end"))      
      
  if (strandSpecific == 1 | strandSpecific == 2) {
        
    if (typeExt == "5") {
      
      # strand +: (1) for each pair of reads, we keep the leftmost 5' end (= 5' cDNA on the + strand)  
      #           (2) we calculate the coverage  
      cov_plus  <- GenomicAlignments::coverage(c(first_ext[BiocGenerics::start(first_ext) <= BiocGenerics::start(last_ext)  & BiocGenerics::strand(first)=="+",], 
                                                 last_ext[BiocGenerics::start(last_ext)   <  BiocGenerics::start(first_ext) & BiocGenerics::strand(first)=="+",]))
      cov_plus  <- cov_plus[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)]  
        
      # strand -: (1) for each pair of reads, we keep the rightmost 5' end (= 5' cDNA on the - strand)  
      #           (2) we calculate the coverage 
      cov_minus <- GenomicAlignments::coverage(c(first_ext[BiocGenerics::start(first_ext) >= BiocGenerics::start(last_ext)  & BiocGenerics::strand(first)=="-",], 
                                                 last_ext[BiocGenerics::start(last_ext)   >  BiocGenerics::start(first_ext) & BiocGenerics::strand(first)=="-",])) 
      cov_minus <- cov_minus[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)]                      
    
    } else { # I assume it is the 3' end
      
      # strand +: (1) for each pair of reads, we keep the rightmost 3' end (= 3' cDNA on the + strand)  
      #           (2) we calculate the coverage  
      cov_plus  <- GenomicAlignments::coverage(c(first_ext[BiocGenerics::start(first_ext) >= BiocGenerics::start(last_ext)  & BiocGenerics::strand(first)=="+",], 
                                                 last_ext[BiocGenerics::start(last_ext)   >  BiocGenerics::start(first_ext) & BiocGenerics::strand(first)=="+",]))
      cov_plus  <- cov_plus[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)]  
        
      # strand -: (1) for each pair of reads, we keep the leftmost 3' end (= 3' cDNA on the - strand)  
      #           (2) we calculate the coverage   
      cov_minus <- GenomicAlignments::coverage(c(first_ext[BiocGenerics::start(first_ext) <= BiocGenerics::start(last_ext)  & BiocGenerics::strand(first)=="-",], 
                                                 last_ext[BiocGenerics::start(last_ext)   <  BiocGenerics::start(first_ext) & BiocGenerics::strand(first)=="-",])) 
      cov_minus <- cov_minus[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)]                          
    }
    
    list(plus = cov_plus, minus = cov_minus)
        
  } else { # I assume it is a single-read library
        
    if (typeExt == "5") {
      
      # by default, we keep the leftmost 5' end and calculate the coverage
      cov_all <- GenomicAlignments::coverage(c(first_ext[BiocGenerics::start(first_ext) <= BiocGenerics::start(last_ext),], 
                                               last_ext[BiocGenerics::start(last_ext)   <  BiocGenerics::start(first_ext),]))
    
    } else { # I assume it is the 3' end
      
      # by default, we keep the rightmost 3' end and calculate the coverage
      cov_all <- GenomicAlignments::coverage(c(first_ext[BiocGenerics::start(first_ext) >= BiocGenerics::start(last_ext),], 
                                               last_ext[BiocGenerics::start(last_ext)   >  BiocGenerics::start(first_ext),]))  
    }
      
    list(all = cov_all[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)])
  }
}

singleReadCov <- function(locus, reads, typeExt, strandSpecific) {
  if (strandSpecific == 2) {
    reads <- BiocGenerics::invertStrand(reads)
  }
      
  ext <- GenomicRanges::resize(GenomicRanges::GRanges(reads), width=1, fix=ifelse(typeExt=="5", "start", "end"))
      
  if (strandSpecific == 1 | strandSpecific == 2) {
        
    cov_by_strand <- lapply(c("+","-"), function(s) { 
      cov <- GenomicAlignments::coverage(ext[BiocGenerics::strand(ext)==s,])
      cov <- cov[[as.character(GenomeInfoDb::seqnames(locus)[1])]][BiocGenerics::start(locus):BiocGenerics::end(locus)]
      cov
    })
  
    names(cov_by_strand) <- c("plus", "minus")
    cov_by_strand
      
  } else {
    list(all = GenomicAlignments::coverage(reads))
  }
}