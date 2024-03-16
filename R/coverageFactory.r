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
  stranded,
  strandSpecific,
  isPairedEnd,
  subsettingBams = FALSE,
  tmpDirectory   = NULL,
  nbThreads      = 1,
  featureCountsOtherParams = list()) {
    
  ## These arguments should not be altered or overridden using 
  ## featureCountsOtherParams list.
  featureCountsOtherParams$files          <- NULL
  featureCountsOtherParams$annot.ext      <- NULL
  featureCountsOtherParams$strandSpecific <- NULL 
  featureCountsOtherParams$nthreads       <- NULL
  featureCountsOtherParams$read2pos       <- NULL
  featureCountsOtherParams$isPairedEnd    <- NULL
  
  ## See oneBinPerBase.
  ## 
  ## MEMORY COMPLEXITY: O(n), n is the locus length 
  annot <- oneBinPerBase(
    locus    = locus,
    stranded = stranded
  )

  ## FeatureCounts goes through all the alignments in a bam to calculate the coverage 
  ## profile of a locus. If this locus is small, it is better to
  ## efficiently build a bam from the subset alignments overlapping this locus (using
  ## samtools view command which directly relies on the index) to then construct the 
  ## associated locus coverage profile. 

  if (subsettingBams) {
    bams  <- createSubBams(
      locus        = locus,
      bams         = bams,
      tmpDirectory = tmpDirectory,
      isPairedEnd  = isPairedEnd
    )
  } else {
    bams  <- list(
      bams  = bams, 
      empty = rep(FALSE, length(bams)) ## (no subsetting) we trust the user
    )
  }

  if (any(!bams$empty)) { ## Run featureCounts if at least one bam is covered.

    featureCounts_args <- c(
      list(
        files             = bams$bams[!bams$empty],
        allowMultiOverlap = TRUE,  
        annot.ext         = annot,
        strandSpecific    = strandSpecific[!bams$empty],
        nthreads          = nbThreads,
        read2pos          = 5, # based on 5' ends of reads
        isPairedEnd       = isPairedEnd[!bams$empty]
      ),
      featureCountsOtherParams
    )

    ## See Rsubread::featureCounts. 
    ##
    ## Rsubread::featureCounts returns, among other outputs, a table of counts 
    ## for each base in each sample.
    ##
    ## PEAK IN MEMORY COMPLEXITY HERE: O(n*s), n is the locus length and s is the 
    ## number of samples 
    five_prime_cov <- do.call(
      Rsubread::featureCounts,
      featureCounts_args
    )

    if (subsettingBams) {
      
      ## Set zero count at each position for samples with empty subset bam.
      if(any(bams$empty)) {
        five_prime_cov <- addZeroCounts(
          bams      = bams$bams,
          empty     = bams$empty,
          fcResults = five_prime_cov
        )
      } 

      ## Remove all temporary subset bams.
      delSubBams(bams = bams$bams, empty = bams$empty)
    }

    ## See formatCoverage. 
    ## 
    ## The counts are compressed using the RLE format. 
    ##
    ## MEMORY COMPLEXITY : We expect to be less than n*s in memory complexity 
    ## after this step since signal is sparse, BUT keep in mind that we are currently 
    ## encountering a peak in memory usage before reaching this point.
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = five_prime_cov
    )
  } else {
    zeroCounts(
      locus    = locus,
      bams     = bams,
      stranded = stranded
    )
  }
}

coverageFactory.threePrime <- function(
  locus,
  bams,
  stranded,
  strandSpecific,
  isPairedEnd,
  subsettingBams = FALSE,
  tmpDirectory   = NULL,
  nbThreads      = 1,
  featureCountsOtherParams = list()) {
    
  ## These arguments should not be altered or overridden using 
  ## featureCountsOtherParams list.
  featureCountsOtherParams$files          <- NULL
  featureCountsOtherParams$annot.ext      <- NULL
  featureCountsOtherParams$strandSpecific <- NULL 
  featureCountsOtherParams$nthreads       <- NULL
  featureCountsOtherParams$read2pos       <- NULL
  featureCountsOtherParams$isPairedEnd    <- NULL
  
  ## See oneBinPerBase.
  ## 
  ## MEMORY COMPLEXITY: O(n), n is the locus length 
  annot <- oneBinPerBase(
    locus    = locus,
    stranded = stranded
  )

  ## FeatureCounts goes through all the alignments in a bam to calculate the coverage 
  ## profile of a locus. If this locus is small, it is better to
  ## efficiently build a bam from the subset alignments overlapping this locus (using
  ## samtools view command which directly relies on the index) to then construct the 
  ## associated locus coverage profile. Only available for unix user with 
  ## installed smatools program.

  if (subsettingBams) {
    bams  <- createSubBams(
      locus        = locus,
      bams         = bams,
      tmpDirectory = tmpDirectory,
      isPairedEnd  = isPairedEnd
    )
  } else {
    bams  <- list(
      bams  = bams, 
      empty = rep(FALSE, length(bams)) ## (no subsetting) we trust the user
    )
  }
  
  if (any(!bams$empty)) { ## Run featureCounts if at least one bam is covered.

    featureCounts_args <- c(
      list(
        files             = bams$bams[!bams$empty],
        allowMultiOverlap = TRUE,  
        annot.ext         = annot,
        strandSpecific    = strandSpecific[!bams$empty],
        nthreads          = nbThreads,
        read2pos          = 3, # based on 3' ends of reads
        isPairedEnd       = isPairedEnd[!bams$empty]
      ),
      featureCountsOtherParams
    )

    ## See Rsubread::featureCounts. 
    ##
    ## Rsubread::featureCounts returns, among other outputs, a table of counts 
    ## for each base in each sample.
    ##
    ## PEAK IN MEMORY COMPLEXITY HERE: O(n*s), n is the locus length and s is the 
    ## number of samples 
    three_prime_cov <- do.call(
      Rsubread::featureCounts,
      featureCounts_args
    )

    if (subsettingBams) {
    
      ## Set zero count at each position for samples with empty subset bam.
      if(any(bams$empty)) {
        three_prime_cov <- addZeroCounts(
          bams      = bams$bams,
          empty     = bams$empty,
          fcResults = three_prime_cov
        )
      } 

      ## Remove all temporary subset bams.
      delSubBams(bams = bams$bams, empty = bams$empty)
    }

    ## See formatCoverage. 
    ## 
    ## The counts are compressed using the RLE format. 
    ##
    ## MEMORY COMPLEXITY : We expect to be less than n*s in memory complexity 
    ## after this step since signal is sparse, BUT keep in mind that we are currently 
    ## encountering a peak in memory usage before reaching this point.
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = three_prime_cov
    )
  } else {
    zeroCounts(
      locus    = locus,
      bams     = bams,
      stranded = stranded
    )
  }
}

coverageFactory.fullLength <- function(
  locus,
  bams,
  stranded,
  strandSpecific,
  isPairedEnd,
  subsettingBams = FALSE,
  tmpDirectory   = NULL,
  nbThreads      = 1,
  featureCountsOtherParams = list()) {
    
  ## These arguments should not be altered or overridden using 
  ## featureCountsOtherParams list.
  featureCountsOtherParams$files          <- NULL
  featureCountsOtherParams$annot.ext      <- NULL
  featureCountsOtherParams$strandSpecific <- NULL 
  featureCountsOtherParams$nthreads       <- NULL
  featureCountsOtherParams$read2pos       <- NULL
  featureCountsOtherParams$isPairedEnd    <- NULL
  
  ## See oneBinPerBase.
  ## 
  ## MEMORY COMPLEXITY: O(n), n is the locus length 
  annot <- oneBinPerBase(
    locus    = locus,
    stranded = stranded
  )

  ## FeatureCounts goes through all the alignments in a bam to calculate the coverage 
  ## profile of a locus. If this locus is small, it is better to
  ## efficiently build a bam from the subset alignments overlapping this locus (using
  ## samtools view command which directly relies on the index) to then construct the 
  ## associated locus coverage profile. Only available for unix user with 
  ## installed smatools program.

  if (subsettingBams) {
    bams  <- createSubBams(
      locus        = locus,
      bams         = bams,
      tmpDirectory = tmpDirectory,
      isPairedEnd  = isPairedEnd
    )
  } else {
    bams  <- list(
      bams  = bams, 
      empty = rep(FALSE, length(bams)) ## (no subsetting) we trust the user
    )
  }
  
  if (any(!bams$empty)) { ## Run featureCounts if at least one bam is covered.
 
    featureCounts_args <- c(
      list(
        files             = bams$bams[!bams$empty],
        allowMultiOverlap = TRUE,  
        annot.ext         = annot,
        strandSpecific    = strandSpecific[!bams$empty],
        nthreads          = nbThreads,
        isPairedEnd       = isPairedEnd[!bams$empty]
      ),
      featureCountsOtherParams
    )

    ## See Rsubread::featureCounts. 
    ##
    ## Rsubread::featureCounts returns, among other outputs, a table of counts 
    ## for each base in each sample.
    ##
    ## PEAK IN MEMORY COMPLEXITY HERE: O(n*s), n is the locus length and s is the 
    ## number of samples 
    fullLength_cov <- do.call(
      Rsubread::featureCounts,
      featureCounts_args
    )

    if (subsettingBams) {
    ## Set zero count at each position for samples with empty subset bam.
      if(any(bams$empty)) {
        fullLength_cov <- addZeroCounts(
          bams      = bams$bams,
          empty     = bams$empty,
          fcResults = fullLength_cov
        )
      } 

      ## Remove all temporary subset bams.
      delSubBams(bams = bams$bams, empty = bams$empty)
    }

    ## See formatCoverage. 
    ## 
    ## The counts are compressed using the RLE format. 
    ##
    ## MEMORY COMPLEXITY : We expect to be less than n*s in memory complexity 
    ## after this step since signal is sparse, BUT keep in mind that we are currently 
    ## encountering a peak in memory usage before reaching this point.  
    formatCoverage(
      locus     = locus,
      stranded  = stranded,
      fcResults = fullLength_cov
    )
  } else {
    zeroCounts(
      locus    = locus,
      bams     = bams,
      stranded = stranded
    )    
  }
}

coverageFactory.average <- function(
  locus,
  bams,
  stranded,
  strandSpecific,
  isPairedEnd,
  subsettingBams = FALSE,
  tmpDirectory   = NULL,
  nbThreads      = 1,
  featureCountsOtherParams = list()) {
    
  ## These arguments should not be altered or overridden using 
  ## featureCountsOtherParams list.
  featureCountsOtherParams$files          <- NULL
  featureCountsOtherParams$annot.ext      <- NULL
  featureCountsOtherParams$strandSpecific <- NULL 
  featureCountsOtherParams$nthreads       <- NULL
  featureCountsOtherParams$read2pos       <- NULL
  featureCountsOtherParams$isPairedEnd    <- NULL
  
  ## See oneBinPerBase.
  ## 
  ## MEMORY COMPLEXITY: O(n), n is the locus length 
  annot <- oneBinPerBase(
    locus    = locus,
    stranded = stranded
  )

  ## FeatureCounts goes through all the alignments in a bam to calculate the coverage 
  ## profile of a locus. If this locus is small, it is better to
  ## efficiently build a bam from the subset alignments overlapping this locus (using
  ## samtools view command which directly relies on the index) to then construct the 
  ## associated locus coverage profile. Only available for unix user with 
  ## installed smatools program.

  if (subsettingBams) {
    bams  <- createSubBams(
      locus        = locus,
      bams         = bams,
      tmpDirectory = tmpDirectory,
      isPairedEnd  = isPairedEnd
    )
  } else {
    bams  <- list(
      bams  = bams, 
      empty = rep(FALSE, length(bams)) ## (no subsetting) we trust the user
    )
  }

  if (any(!bams$empty)) { ## Run featureCounts if at least one bam is covered.

    ends_cov <- lapply(c(3,5), function(end) {
      featureCounts_args <- c(
        list(
          files             = bams$bams[!bams$empty],
          allowMultiOverlap = TRUE,  
          annot.ext         = annot,
          strandSpecific    = strandSpecific[!bams$empty],
          nthreads          = nbThreads,
          read2pos          = end, # based on 3' or 5' ends of reads
          isPairedEnd       = isPairedEnd[!bams$empty]
        ),
        featureCountsOtherParams
      )

      ## See Rsubread::featureCounts. 
      ##
      ## Rsubread::featureCounts returns, among other outputs, a table of counts 
      ## for each base in each sample.
      ##
      ## PEAK IN MEMORY COMPLEXITY HERE: O(n*s), n is the locus length and s is the 
      ## number of samples 
      one_end_cov <- do.call(
        Rsubread::featureCounts,
        featureCounts_args
      )

      if (subsettingBams) {
      ## Set zero count at each position for samples with empty subset bam.
        if(any(bams$empty)) {
          one_end_cov <- addZeroCounts(
            bams      = bams$bams,
            empty     = bams$empty,
            fcResults = one_end_cov
          )
        } 
      }

      ## See formatCoverage. 
      ## 
      ## The counts are compressed using the RLE format. 
      ##
      ## MEMORY COMPLEXITY : We expect to be less than n*s in memory complexity 
      ## after this step since signal is sparse, BUT keep in mind that we are currently 
      ## encountering a peak in memory usage before reaching this point.    
      formatCoverage(
        locus     = locus,
        stranded  = stranded,
        fcResults = one_end_cov
      )
    })

    if (subsettingBams) {   
      ## Remove all temporary subset bams.
      delSubBams(bams = bams$bams, empty = bams$empty)
    }

    names(ends_cov) <- c("3","5")

    ## calculate the geometric mean of 3' and 5' coverages
    coverage <- list()
    for (sample in names(ends_cov[["3"]])) {
      for (strand in names(ends_cov[["3"]][[sample]])) {
        coverage[[sample]][[strand]] <- exp((log(ends_cov[["3"]][[sample]][[strand]]+1)+
          log(ends_cov[["5"]][[sample]][[strand]]+1))/2)-1
      }
    }
    
    coverage
  } else {
    zeroCounts(
      locus    = locus,
      bams     = bams,
      stranded = stranded
    )
  }
}

formatCoverage <- function(
  locus, 
  stranded,
  fcResults) {

  ## formatCoverage generates a structured list of compressed coverages 
  ## (RLE format). Each entry in this list represents a specific 'sample', 
  ## the coverage data is further organized by 'strand', referring to the two
  ## strands.

  coverages <- list()
  relative_start <- 1
  relative_end   <- GenomicRanges::width(locus)
  
  if (stranded) {
    for (i in seq_along(fcResults$targets)) {

      ## coverage for both plus and minus strands is sequentially organized in
      ## the table returned by featureCounts
      coverages[[fcResults$target[[i]]]][["plus"]] <- S4Vectors::Rle(fcResults$counts[
        relative_start:relative_end,i
      ])
      
      coverages[[fcResults$target[[i]]]][["minus"]] <- S4Vectors::Rle(fcResults$counts[
        (relative_end+1):nrow(fcResults$counts),i
      ])
    }
  } else {
    for (i in seq_along(fcResults$targets)) {
      coverages[[fcResults$target[[i]]]][["all"]] <- S4Vectors::Rle(fcResults$counts[
        relative_start:relative_end,i
      ])
    }
  }

  coverages
}

oneBinPerBase <- function(
  locus, 
  stranded) {

  ## oneBinPerBase generates an annotation for each base within the locus specified by
  ## the user. These annotations are systematically stored in a 'data.frame',
  ## adhering to the Simple Annotation Format (SAF). The data.frame includes
  ## the following columns: 'GeneID', 'Chr', 'Start', 'End', and 'Strand'.

  if (stranded) {
    GeneID <- c( 
      paste0("+" , GenomicRanges::start(locus):GenomicRanges::end(locus)),
      paste0("-" , GenomicRanges::start(locus):GenomicRanges::end(locus))
    )
    Strand <- rep(c("+" ,"-"), each=GenomicRanges::width(locus))
  } else {

    ## In cases where the reads are unstranded, they are counted onto the 
    ## forward strand by default. Does it even work ?
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


createSubBams <- function(
  locus, 
  bams,
  isPairedEnd,
  tmpDirectory) { 

  ## createSubBams generates subset BAMs containing only alignments overlapping 
  ## the current locus. It calls the samtools view program.
  sub_bams <- list(bams=c(), empty=c())
  for (i_bam in seq_along(bams)) {
    
    bam <- bams[[i_bam]]

    param <- Rsamtools::ScanBamParam(
      which = locus, 
      what  = c("rname", "strand", "pos", "qwidth", "seq")
    )

    if (isPairedEnd[[i_bam]]) {
      suppressWarnings(alignments <- GenomicAlignments::readGAlignmentPairs(
        file  = bam, 
        param = param
      ))
    } else {
      suppressWarnings(alignments <- GenomicAlignments::readGAlignments(
        file  = bam, 
        param = param
      ))
    }

    sub_bam <- paste0(
      as.character(GenomeInfoDb::seqnames(locus)),
      "_",
      GenomicRanges::start(locus),
      "_",
      GenomicRanges::end(locus),
      "_",
      basename(bam)
    )
    
    sub_bam <- file.path(tmpDirectory, sub_bam)

    sub_bams$bams <- c(sub_bams$bams, sub_bam)

    if (length(alignments)>0) {
      sub_bams$empty <-  c(sub_bams$empty, FALSE)
      rtracklayer::export(alignments , Rsamtools::BamFile(sub_bam))
    } else {
      sub_bams$empty <-  c(sub_bams$empty, TRUE)
    }
  }

  sub_bams
}

delSubBams <- function(bams, empty) {

  ## delSubBams removes temporary subset BAMs.

  unlink(bams[!empty])
  unlink(paste0(bams[!empty], ".bai"))
}

addZeroCounts <- function(
  bams,
  empty,
  fcResults) {

  ## addZeroCounts sets zero count at each position for samples with empty 
  ## subset bam.

  tmp <- as.data.frame(lapply(
    basename(bams[empty]), 
    function(x) rep(0, nrow(fcResults$counts))
  ))

  names(tmp)        <- basename(bams[empty])
  row.names(tmp)    <- row.names(fcResults$counts)
  fcResults$counts  <- cbind(fcResults$counts, tmp)
  fcResults$counts  <- fcResults$counts[,basename(bams)]
  fcResults$targets <- basename(bams)
  fcResults
}


zeroCounts <- function(
  locus, 
  bams, 
  stranded) {

  cov <- lapply(bams$bams, function(bam) {
    if (stranded) {
      list(
        plus  = S4Vectors::Rle(rep(0, GenomicRanges::width(locus))),
        minus = S4Vectors::Rle(rep(0, GenomicRanges::width(locus)))
      )  
    } else {
      list(all = S4Vectors::Rle(rep(0, GenomicRanges::width(locus))))
    }
  })
  names(cov) <- basename(bams$bams)
  cov
}