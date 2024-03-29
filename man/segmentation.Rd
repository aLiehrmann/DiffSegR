% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation.r
\name{segmentation}
\alias{segmentation}
\title{Summarize the differential transcription landscape}
\usage{
segmentation(
  data,
  weightType = "unweighted",
  modelSelectionType = "yao",
  featureCountsType = "fromBam",
  compressed = TRUE,
  alpha = 2,
  segmentNeighborhood = FALSE,
  Kmax = NULL,
  verbose = FALSE,
  nbThreadsGridSearch = 1,
  alphas = NULL,
  gridSearch = FALSE,
  outputDirectory = ".",
  nbThreadsFeatureCounts = 1,
  strandSpecific = 1,
  read2pos = NULL,
  isPairedEnd = FALSE,
  featureCountsOtherParams = list()
)
}
\arguments{
\item{data}{The \code{List} object returned by \code{\link[=loadData]{loadData()}}.}

\item{weightType}{A \code{String}. Select the type of weights associated to
the log2-FC per-base:
\itemize{
\item unweighted   : all observations have the same weight (=1) ;
\item zeroInflated : low counts have less weight.
}}

\item{modelSelectionType}{A \code{String}. Select the penalty used by FPOP:
\itemize{
\item yao : Yao's penalty \code{alpha*sigma^2*log(n)}.
}}

\item{featureCountsType}{A \code{String}. Select how to summarize counts:
\itemize{
\item fromBam      : from bam files using \link[Rsubread:featureCounts]{Rsubread::featureCounts} ;
\item fromCoverage : from coverage profiles.
}}

\item{compressed}{A \code{Boolean}. Indicate if the observations have to be
compressed (does not change the segmentation results and decreases
the running time).}

\item{alpha}{A \code{Double}. Segmentation hyperparameter used in Yao's penalty:
\code{alpha*sigma^2*log(n)}. The number of changepoints returned by
\code{\link[fpopw:Fpop_w]{fpopw::Fpop_w()}} is a decreasing function of \code{alpha}.}

\item{segmentNeighborhood}{A \code{Boolean}. Indicate the weighted pDPA algorithm
(Rigaill 2010 and 2015) has to be used to explore the space of segmentations.}

\item{Kmax}{An \code{Integer}. Segmentations with 1 to Kmax segments are recovered
using segment neighborhood.}

\item{verbose}{A \code{Boolean}. Should all the operations performed be
displayed ?}

\item{nbThreadsGridSearch}{An \code{Integer}. Number of threads used by the grid
search procedure.}

\item{alphas}{A vector of \code{Double}. A series of alphas used by the grid
search procedure.}

\item{gridSearch}{A \code{Boolean}. Indicate if a grid search has to be used
to explore the space of segmentations.}

\item{outputDirectory}{A \code{String}. Absolute path to the output directory.}

\item{nbThreadsFeatureCounts}{An \code{Integer}. Passed on to
\code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.}

\item{strandSpecific}{An \code{Integer}. Passed on to \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.}

\item{read2pos}{An \code{Integer}. Passed on to \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.}

\item{isPairedEnd}{Passed on to \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.}

\item{featureCountsOtherParams}{A \code{List}. Other paramters passed on to
\code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.}
}
\value{
A \code{SummarizedExperiment} object that contains region boundaries and
associated count matrix.
}
\description{
In brief, \code{segmentation()} (i) identifies segments with homogeneous log2-FC
along the locus of interest and (ii) counts reads overlapping them.

In details, \code{segmentation()} (i) starts by calling \code{\link[fpopw:Fpop_w]{fpopw::Fpop_w()}} on the
log2-FC profile of both strands to retreive the segments boundaries
(changepoints). The hyperparameter \code{alpha} specified by the user controls
the number of returned segments. Corresponding genomic regions are stored as
a GenomicRanges object. (ii) Finally, \code{segmentation()} calls
\code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}} which assigns to the regions the mapped reads
from each replicate of each biological condition. By default a read is
allowed to be assigned to more than one region if it is found to overlap
with more than one region. \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}} ends by building a
count matrix. \code{segmentation()} ends by returning the regions and associated
count matrix as \code{SummarizedExperiment} object.
}
\examples{

# Create a working directory for running the example.
working_directory <- "./DIFFSEGR_TEST"
dir.create(working_directory)

# Save sample information in a text file.
sample_info <- data.frame(
  sample    = c("pnp1_1_1", "pnp1_1_2", "wt_1", "wt_2"),
  condition = rep(c("pnp1_1", "wt"), each = 2),
  replicate = rep(1:2,2),
  bam       = sapply(
    c("pnp1_1_1_ChrC_71950_78500.bam", 
      "pnp1_1_2_ChrC_71950_78500.bam",
      "wt_1_ChrC_71950_78500.bam",
      "wt_2_ChrC_71950_78500.bam"
     ),
     function(bam) system.file("extdata", bam, package = "DiffSegR")
  ),
  coverage  = file.path(
    working_directory,
    paste0(c("pnp1_1_1", "pnp1_1_2", "wt_1", "wt_2"), ".rds")
  )
)
write.table(
  sample_info, 
  file.path(working_directory, "sample_info.txt")
)

# Build coverages and log2-FC.
data <- loadData(
  sampleInfo         = file.path(working_directory,"sample_info.txt"),
  locus        = list(
    seqid      = "ChrC", 
    chromStart = 71950, 
    chromEnd   = 78500
  ),
  referenceCondition = "wt",
  stranded           = TRUE,
  fromBam            = TRUE,
  nbThreads          = 1
)

# Summarize the differential landscape.
SExp <- segmentation(
  data                   = data, 
  nbThreadsFeatureCounts = 1,
  outputDirectory        = working_directory
)

# In genomic order, first to fifth regions and ...
print(SummarizedExperiment::mcols(SExp)[1:5,])
# ... associated counts.
print(SummarizedExperiment::assay(SExp)[1:5,])

# delete working directory 
unlink(working_directory, recursive = TRUE)

}
