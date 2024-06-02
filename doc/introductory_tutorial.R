## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
working_directory <- "./DIFFSEGR_TEST"
dir.create(working_directory, showWarnings = FALSE)

## -----------------------------------------------------------------------------
#- create sample information table --------------------------------------------#
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

#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

## -----------------------------------------------------------------------------
library(DiffSegR)
nb_threads_tot = 4
nb_threads_locus = 4

## -----------------------------------------------------------------------------
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

print(data)

## -----------------------------------------------------------------------------
coverage(
  data         = data,
  coverageType = "average"
)

## ----eval=FALSE---------------------------------------------------------------
#  coverage(
#    data           = data,
#    coverageType   = "average",
#    subsettingBams = TRUE
#  )

## -----------------------------------------------------------------------------
features <- segmentationLFC(
  data  = data, 
  alpha = 2
)

## -----------------------------------------------------------------------------
SExp <- counting(
  data     = data, 
  features = features
)

## -----------------------------------------------------------------------------
#- number of segments found ---------------------------------------------------# 
SummarizedExperiment::strand(SExp)
#- display first five segments ------------------------------------------------#
segment_coordinates <- as.data.frame(SummarizedExperiment::rowRanges(SExp))
knitr::kable(segment_coordinates[1:5, colnames(segment_coordinates)%in%c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand")
]) 
#- display counts associated to first five segments ---------------------------#
counts <- SummarizedExperiment::assay(SExp)
knitr::kable(counts[1:5,])

## -----------------------------------------------------------------------------
dds <- dea(
  SExp              = SExp,
  design            = ~ condition,
  sizeFactors       = rep(1,4),
  significanceLevel = 1e-2
)

## -----------------------------------------------------------------------------
#- number of DERs -------------------------------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$DER,]
length(DERs)
#- display first five DERs by position ----------------------------------------#
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))
knitr::kable(DERs[1:5, colnames(DERs)%in%c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand",
  "log2FoldChange",
  "padj")
])

## -----------------------------------------------------------------------------
exons_path <- system.file(
 "extdata", 
 "TAIR10_ChrC_exon.gff3", 
 package = "DiffSegR"
)
annotated_DERs <- annotateNearest(
  features = SummarizedExperiment::rowRanges(dds),   
  annotations = exons_path,
  select      = c("seqnames", "start", "end", "width", "strand", 
     "description", "distance", "Parent", "log2FoldChange", "pvalue"),
  orderBy     = "pvalue",
  outputDirectory = working_directory,
)
knitr::kable(annotated_DERs, row.names = FALSE, digits = 3)

## -----------------------------------------------------------------------------
ChrC_path  <- system.file(
 "extdata", 
 "TAIR10_ChrC.fa", 
 package = "DiffSegR"
)
exportResults(
 data            = data,
 features        = SummarizedExperiment::rowRanges(dds),
 sizeFactors     = DESeq2::sizeFactors(dds),
 outputDirectory = working_directory,
 genome          = ChrC_path,
 annotation      = exons_path,
 select          = c("log2FoldChange", "pvalue", "padj")
) 

## ----echo=FALSE, out.width='100%', fig.cap="DiffSegR output between positions 76,990 to 78,248 on the chloroplast genome. The session is loaded in IGV 2.17.2 for Linux."----
knitr::include_graphics('igv_snapshot.png')

## ----eval=FALSE---------------------------------------------------------------
#  setPaths(
#    directory = "collaborator_igv_session(s)_directory",
#    genome    = "collaborator_genome_path",
#    annotations = "collaborator_annotations_path"
#  )

