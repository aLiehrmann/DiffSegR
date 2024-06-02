## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----include = FALSE----------------------------------------------------------
working_directory <- "./DIFFSEGR_TEST"
dir.create(working_directory, showWarnings = FALSE)
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
  isPairedEnd    = rep(FALSE, 4),
  strandSpecific = rep(1, 4)
)
library(DiffSegR)
nb_threads_tot = 4
nb_threads_locus = 4
data <- newExperiment(
  sampleInfo   = sample_info,
  loci         = data.frame(
    seqid      = "ChrC", 
    chromStart = 71950, 
    chromEnd   = 78500,
    locusID    = "psbB_petD"
  ),
  referenceCondition = "wt",
  otherCondition     = "pnp1_1",
  nbThreads          = nb_threads_tot,
  nbThreadsByLocus   = nb_threads_locus,
  coverage           = working_directory
)
coverage(
  data         = data,
  coverageType = "average"
)

## -----------------------------------------------------------------------------
features <- segmentationLFC(
  data                = data, 
  outputDirectory     = working_directory,
  segmentNeighborhood = TRUE,
  Kmax                = 50
)

## ----echo=FALSE, dpi=300, out.width='100%', fig.align='center', fig.cap="Path of all reachable segmentations up to 50 segments of the log2-FC on reverse strand. Corresponding alpha values on the log scale can be read on X-axis."----
all_models <- readRDS("DIFFSEGR_TEST/psbB_petD/all_models_minus.rds")
ggplot2::ggplot(
  all_models,
  ggplot2::aes(
    x = log2(alphas), 
    y = K
  )
)+
ggplot2::geom_line() +
ggplot2::geom_point(size=1) +
ggplot2::geom_vline(
  xintercept = log2(5.46), 
  color      = "red"
)+
ggplot2::geom_hline(
  yintercept = 11, 
  color      = "blue"
)+
ggplot2::xlab("log2(alpha)")+
ggplot2::ylab("number of segments (K)")+
ggplot2::theme(
  text = ggplot2::element_text(size = 20)
) +
ggplot2::theme_bw()

## ----include = FALSE----------------------------------------------------------
SExp <- counting(
  data     = data, 
  features = features
)

## -----------------------------------------------------------------------------
dds <- dea(
  SExp                      = SExp,
  design                    = ~ condition,
  sizeFactors               = rep(1,4),
  predicate                 = function(x) 2^abs(x$log2FoldChange)>2,
  postHoc_significanceLevel = 0.05,
  postHoc_tdpLowerBound     = 0.95
)

#- check fold-change of differentially expressed segments ---------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$DER,]
2^abs(SummarizedExperiment::mcols(DERs)$log2FoldChange)

