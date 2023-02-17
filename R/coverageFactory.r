coverageFactory <- function(type) {
  fn <- switch(type,
    fivePrime   = coverageFactory.fivePrime,
    threePrime  = coverageFactory.threePrime,
    average     = coverageFactory.average,
    fullLength  = coverageFactory.fullLength,
    center  = coverageFactory.center,
    eval(parse(text=type))
  )
  attr(fn, "type") <- type
  fn
}

coverageFactory.fivePrime <- function(
  locus,
  stranded,
  bams,
  strandSpecific = 1,
  nbThreads      = 1,
  readLength     = 1,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
    
    featureCountsOtherParams$files          <- NULL
    featureCountsOtherParams$annot.ext      <- NULL
    featureCountsOtherParams$strandSpecific <- NULL 
    featureCountsOtherParams$nthreads       <- NULL
    featureCountsOtherParams$read2pos       <- NULL
    featureCountsOtherParams$isPairedEnd    <- NULL
    
    if (!stranded) {
      strandSpecific <- 0
    }
    
    annot <- oneBinPerBase(
      locus = locus,
      stranded = stranded
    )
    
    five_prime_ends <- do.call(
      Rsubread::featureCounts,
      c(
        list(
          files          = bams,
          ##- featureCounts cuts reads before counting -----------------------##
          #allowMultiOverlap = TRUE,  
          annot.ext      = annot,
          strandSpecific = strandSpecific,
          nthreads       = nbThreads,
          read2pos       = 5,
          isPairedEnd    = isPairedEnd
        ),
        featureCountsOtherParams
      )
    )
    
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = five_prime_ends
    )
}

coverageFactory.threePrime <- function(
  locus,
  stranded,
  bams,
  strandSpecific = 1,
  nbThreads      = 1,
  readLength     = 1,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
    
    featureCountsOtherParams$files          <- NULL
    featureCountsOtherParams$annot.ext      <- NULL
    featureCountsOtherParams$strandSpecific <- NULL 
    featureCountsOtherParams$nthreads       <- NULL
    featureCountsOtherParams$read2pos       <- NULL
    featureCountsOtherParams$isPairedEnd    <- NULL
    
    if (!stranded) {
      strandSpecific <- 0
    }
    
    annot <- oneBinPerBase(
      locus    = locus,
      stranded = stranded
    )
    
    three_prime_ends <- do.call(
      Rsubread::featureCounts,
      c(
        list(
          files          = bams,
          ##- featureCounts cuts reads before counting -----------------------## 
          #allowMultiOverlap = TRUE,  
          annot.ext      = annot,
          strandSpecific = strandSpecific,
          nthreads       = nbThreads,
          read2pos       = 3,
          isPairedEnd    = isPairedEnd
        ),
        featureCountsOtherParams
      )
    )
    
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = three_prime_ends
    )
}

coverageFactory.center <- function(
  locus,
  stranded,
  bams,
  strandSpecific = 1,
  nbThreads      = 1,
  readLength     = 1,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
  
    featureCountsOtherParams$files          <- NULL
    featureCountsOtherParams$annot.ext      <- NULL
    featureCountsOtherParams$strandSpecific <- NULL 
    featureCountsOtherParams$nthreads       <- NULL
    featureCountsOtherParams$read2pos       <- NULL
    featureCountsOtherParams$isPairedEnd    <- NULL
    featureCountsOtherParams$readShiftType  <- NULL
    featureCountsOtherParams$readShiftSize  <- NULL
    
    if (!stranded) {
      strandSpecific <- 0
    }
    
    annot <- oneBinPerBase(
      locus    = locus,
      stranded = stranded
    )
    
    center <- do.call(
      Rsubread::featureCounts,
      c(
        list(
          files          = bams,
          ##- featureCounts cuts reads before counting -----------------------## 
          #allowMultiOverlap = TRUE,  
          annot.ext      = annot,
          strandSpecific = strandSpecific,
          nthreads       = nbThreads,
          readShiftType  = "right",
          readShiftSize  = floor(readLength/2),
          read2pos       = 5,
          isPairedEnd    = isPairedEnd
        ),
        featureCountsOtherParams
      )
    ) 
    
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = center
    )
}

coverageFactory.fullLength <- function(
  locus,
  stranded,
  bams,
  strandSpecific = 1,
  nbThreads      = 1,
  readLength     = 1,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
  
    featureCountsOtherParams$files             <- NULL
    featureCountsOtherParams$annot.ext         <- NULL
    featureCountsOtherParams$allowMultiOverlap <- NULL
    featureCountsOtherParams$strandSpecific    <- NULL 
    featureCountsOtherParams$nthreads          <- NULL
    featureCountsOtherParams$read2pos          <- NULL
    featureCountsOtherParams$isPairedEnd       <- NULL
    
    if (!stranded) {
      strandSpecific <- 0
    }
    
    annot <- oneBinPerBase(
      locus    = locus,
      stranded = stranded
    )
    
    full_length <- do.call(
      Rsubread::featureCounts,
      c(
        list(
          files             = bams,
          allowMultiOverlap = TRUE,  
          annot.ext         = annot,
          strandSpecific    = strandSpecific,
          nthreads          = nbThreads,
          read2pos          = NULL,
          isPairedEnd       = isPairedEnd
        ),
        featureCountsOtherParams
      )
    )
    
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = full_length
    )
}

coverageFactory.average <- function(
  locus,
  stranded,
  bams,
  strandSpecific = 1,
  nbThreads      = 1,
  readLength     = 1,
  isPairedEnd    = FALSE,
  featureCountsOtherParams = list()) {
  
    featureCountsOtherParams$files          <- NULL
    featureCountsOtherParams$annot.ext      <- NULL
    featureCountsOtherParams$strandSpecific <- NULL 
    featureCountsOtherParams$nthreads       <- NULL
    featureCountsOtherParams$read2pos       <- NULL
    featureCountsOtherParams$isPairedEnd    <- NULL
    
    if (!stranded) {
      strandSpecific <- 0
    }
    
    annot <- oneBinPerBase(
      locus    = locus,
      stranded = stranded
    )
    
    three_prime_ends <- do.call(
      Rsubread::featureCounts,
      c(
        list(
          files          = bams,
          ##- featureCounts cuts reads before counting -----------------------##
          #allowMultiOverlap = TRUE,  
          annot.ext      = annot,
          strandSpecific = strandSpecific,
          nthreads       = nbThreads,
          read2pos       = 3,
          isPairedEnd    = isPairedEnd
        ),
        featureCountsOtherParams
      )
    )
    
    five_prime_ends <- do.call(
      Rsubread::featureCounts,
      c(
        list(
          files          = bams,
          ##- featureCounts cuts reads before counting -----------------------##
          #allowMultiOverlap = TRUE,  
          annot.ext      = annot,
          strandSpecific = strandSpecific,
          nthreads       = nbThreads,
          read2pos       = 5,
          isPairedEnd    = isPairedEnd
        ),
        featureCountsOtherParams
      )
    )
    
    three_prime_ends_cov <- formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = three_prime_ends
    )
    
    five_prime_ends_cov <- formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = five_prime_ends
    )
    
    ##- compute geometric mean of 3' and 5' ends ------------------------------##
    coverage <- list()
    for (sample in names(three_prime_ends_cov)) {
      for (strand in names(three_prime_ends_cov[[sample]])) {
        cov <- matrix(
          c(
            three_prime_ends_cov[[sample]][[strand]],
            five_prime_ends_cov[[sample]][[strand]]
          ),
          nrow = length(three_prime_ends_cov[[sample]][[strand]])
        )
        coverage[[sample]][[strand]] <- exp(rowMeans(log(cov+1)))-1
        coverage[[sample]][[strand]] <- S4Vectors::Rle(
          coverage[[sample]][[strand]] 
        )
     }
   }
   coverage
}

formatCoverage <- function(
  locus, 
  stranded,
  fcResults) {
  coverages <- list()
  relative_start <- 1
  relative_end   <- GenomicRanges::width(locus)
  if (stranded) {
    for (i in seq_along(fcResults$targets)){
      coverages[[fcResults$target[[i]]]][["plus"]] <- fcResults$counts[
        relative_start:relative_end,i
      ]
      coverages[[fcResults$target[[i]]]][["plus"]] <- S4Vectors::Rle(
        coverages[[fcResults$target[[i]]]][["plus"]]
      )
      
      coverages[[fcResults$target[[i]]]][["minus"]] <- fcResults$counts[
        (relative_end+1):nrow(fcResults$counts),i
      ]
      coverages[[fcResults$target[[i]]]][["minus"]] <- S4Vectors::Rle(
        coverages[[fcResults$target[[i]]]][["minus"]]
      )
    }
  } else {
    for (i in seq_along(fcResults$targets)){
      coverages[[fcResults$target[[i]]]][["all"]] <- fcResults$counts[
        relative_start:relative_end,i
      ]
      coverages[[fcResults$target[[i]]]][["all"]] <- S4Vectors::Rle(
        coverages[[fcResults$target[[i]]]][["all"]]
      )
    }
  }
  coverages
}

oneBinPerBase <- function(locus, stranded) {
  if (stranded) {
    GeneID <- c( 
      paste0("+" , GenomicRanges::start(locus):GenomicRanges::end(locus)),
      paste0("-" , GenomicRanges::start(locus):GenomicRanges::end(locus))
    )
    Strand <- rep(c("+" ,"-"), each=GenomicRanges::width(locus))
  } else {
    GeneID <- GenomicRanges::start(locus):GenomicRanges::end(locus)
    Strand <- "+"
  }
  
  data.frame(
    GeneID = GeneID,
    Chr    = as.character(GenomeInfoDb::seqnames(locus)),
    Start  = GenomicRanges::start(locus):GenomicRanges::end(locus),
    End    = GenomicRanges::start(locus):GenomicRanges::end(locus),
    Strand = Strand
  ) 
}
