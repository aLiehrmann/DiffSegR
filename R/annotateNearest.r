#' Annotate each provided differentialy expressed region (DER) with its
#' closest annotation(s)
#'
#' @description 
#' `annotateNearest()` relies on user specified annotations in gff3 or gtf 
#' format to annotate the DERs found during the DEA. Seven classes of labels 
#' translate the relative positions of the DER and its closest annotation(s): 
#' * antisense ; 
#' * upstream ; 
#' * downstream ; 
#' * inside ; 
#' * overlapping 3’ ; 
#' * overlapping 5’ ; 
#' * overlapping 5’, 3’.
#' 
#' @param data The `List` object returned by [DiffSegR::loadData()].
#' @param dds The `DESeqDataSet` object returned by [DiffSegR::dea()].
#' @param annotations A `String`. Path to the annotation file (in gff or gtf 
#' format).
#' @param annotationsType A vector of `String`. Keep only matching annotations.
#' @param orderBy A `String`. The name of the column on which returned 
#' annotations are sorted.
#' @param select A vector of `String`. Keep only columns matching these 
#' names in the returned `Data.frame`. By default all columns from 
#' mcols(dds) and mcols(annotations) are kept.
#' @returns Annotated differentialy expressed regions as `Data.frame` object.
#' 
#' @examples
#' # Create a working directory for running the example.
#' working_directory <- "./DIFFSEGR_TEST"
#' dir.create(working_directory)
#' 
#' # Save sample information in a text file.
#' sample_info <- data.frame(
#'   sample    = c("pnp1_1_1", "pnp1_1_2", "wt_1", "wt_2"),
#'   condition = rep(c("pnp1_1", "wt"), each = 2),
#'   replicate = rep(1:2,2),
#'   bam       = sapply(
#'     c("pnp1_1_1_ChrC_71950_78500.bam", 
#'       "pnp1_1_2_ChrC_71950_78500.bam",
#'       "wt_1_ChrC_71950_78500.bam",
#'       "wt_2_ChrC_71950_78500.bam"
#'      ),
#'      function(bam) system.file("extdata", bam, package = "DiffSegR")
#'   ),
#'   coverage  = file.path(
#'     working_directory,
#'     paste0(c("pnp1_1_1", "pnp1_1_2", "wt_1", "wt_2"), ".rds")
#'   )
#' )
#' write.table(
#'   sample_info, 
#'   file.path(working_directory, "sample_info.txt")
#' )
#' 
#' # Build coverages and log2-FC per-base.
#' data <- loadData(
#'   sampleInfo         = file.path(working_directory,"sample_info.txt"),
#'   locus        = list(
#'     seqid      = "ChrC", 
#'     chromStart = 71950, 
#'     chromEnd   = 78500
#'   ),
#'   referenceCondition = "wt",
#'   stranded           = TRUE,
#'   fromBam            = TRUE,
#'   nbThreads          = 1
#' )
#' 
#' # Summarize the differential landscape.
#' SExp <- segmentation(
#'   data                   = data, 
#'   nbThreadsFeatureCounts = 1,
#'   outputDirectory        = working_directory
#' )
#' 
#' # Differential expression analysis. We control all the resulting 
#' # p-values using a Benjamini–Hochberg procedure at level 0.01.
#' dds <- dea(
#'   data              = data,
#'   SExp              = SExp,
#'   design            = ~ condition,
#'   sizeFactors       = rep(1,4),
#'   significanceLevel = 1e-2
#' )
#' 
#' # Annotate each differential expressed region with its closest annotation(s).
#' exons_path <- system.file(
#'   "extdata", 
#'   "TAIR10_ChrC_exon.gff3", 
#'   package = "DiffSegR"
#' )
#' annotated_DERs <- annotateNearest(
#'   data        = data,
#'   dds         = dds,
#'   annotations = exons_path,
#'   select      = c("seqnames", "start", "end", "width", "strand", 
#'     "description", "distance", "Parent", "log2FoldChange", "pvalue"),
#'   orderBy     = "pvalue"
#' )
#' print(annotated_DERs[1:10,])
#' 
#' # delete working directory 
#' unlink(working_directory, recursive = TRUE)
#' 
#' @export
annotateNearest <- function(
  data,
  dds, 
  annotations,
  annotationsType = NULL,
  orderBy = "start",
  select  = NULL) {
  
  dds <- dds[GenomicRanges::mcols(dds)$rejectedHypotheses==TRUE,]
  if (length(dds)<1){
    stop("No differential expressed region to annotate.")
  }
  
  ##- loading annotations ----------------------------------------------------##
  message(
    " > loading annotations from ",
    annotations,
    " ..."
  )
  annotations <- rtracklayer::import(annotations)
  if (!is.null(annotationsType)) {
    annotations <- annotations[annotations$type %in% annotationsType,]
  }
  annotations <- IRanges::subsetByOverlaps(annotations,data$locus)
  
  ##- (i) DERs and annotations on same strand --------------------------------##
  message(" > annotate differential expressed regions ...")
  annotated_DERs <- findNearestAnnotation(
    dds         = dds,
    annotations = annotations,
    antisense   = FALSE
  )
  
  ##- (ii) overlapping DERs and annotations on opposite strand (antisense) ---##
  annotated_DERs_as <- findNearestAnnotation(
    dds         = dds,
    annotations = annotations,
    antisense   = TRUE
  )
  ##- merge (i) and (ii) -----------------------------------------------------##
  all_annotated_DERs <- rbind(annotated_DERs, annotated_DERs_as)

  ##- labeling the annotations ---------------------------------##
  all_annotated_DERs <- labeling(
    all_annotated_DERs = all_annotated_DERs
  )
  all_annotated_DERs <- data.frame(lapply(all_annotated_DERs ,unlist))
  ##- keep only selected columns ---------------------------------------------##
  if (is.null(select)) {
    all_annotated_DERs[order(all_annotated_DERs[[orderBy]]),]
  } else {
    all_annotated_DERs[order(all_annotated_DERs[[orderBy]]),select]
  }
}

findNearestAnnotation <- function(
  dds,
  annotations,
  antisense = FALSE) {
  
  association_table <- GenomicRanges::distanceToNearest(
    x             = SummarizedExperiment::rowRanges(dds), 
    subject       = annotations,
    select        = "all",
    ignore.strand = antisense
  )  
  
  nearest_annotations <- as.data.frame(annotations[association_table@to,])
  nearest_annotations$startAnnotation  <- nearest_annotations$start
  nearest_annotations$endAnnotation    <- nearest_annotations$end
  nearest_annotations$strandAnnotation <- nearest_annotations$strand
  nearest_annotations$widthAnnotation  <- nearest_annotations$width
  nearest_annotations <- nearest_annotations[
    ,!names(nearest_annotations)%in%c("start","end","strand","width")
  ]
  
  annotated_DERs <- cbind(
    as.data.frame(
      SummarizedExperiment::rowRanges(dds[association_table@from,]),
      row.names = NULL
    ),
    nearest_annotations,
    distance = GenomicRanges::mcols(association_table)$distance
  )
  
  if (antisense) {
    annotated_DERs <- annotated_DERs[
      annotated_DERs$strand != annotated_DERs$strandAnnotation & 
        annotated_DERs$distance == 0,
    ]
  }
  
  annotated_DERs
}

labeling <- function(
  all_annotated_DERs) {
  do.call(rbind, lapply(
    1:nrow(all_annotated_DERs), 
    function(i){
      nearest_annotation <- all_annotated_DERs[i,]
      if (nearest_annotation$strand != nearest_annotation$strandAnnotation & 
          nearest_annotation$distance == 0){
        nearest_annotation$description <- "antisense"
      } else if (nearest_annotation$end<nearest_annotation$startAnnotation){
        if (nearest_annotation$strand == "+"){
          nearest_annotation$description <- "upstream"
        } else {
          nearest_annotation$description <- "downstream"
        }
      } else if (nearest_annotation$start>nearest_annotation$endAnnotation) {
        if (nearest_annotation$strand == "+"){
          nearest_annotation$description <- "downstream"
        } else {
          nearest_annotation$description <- "upstream"
        }
      } else if (nearest_annotation$start>=nearest_annotation$startAnnotation & 
                 nearest_annotation$end<=nearest_annotation$endAnnotation) {
        nearest_annotation$description <- "inside"
      } else if (nearest_annotation$start<=nearest_annotation$startAnnotation & 
                 nearest_annotation$end>=nearest_annotation$endAnnotation){
        nearest_annotation$description <- "overlapping 5', 3'"
      } else if (nearest_annotation$start<=nearest_annotation$endAnnotation & 
                 nearest_annotation$start>nearest_annotation$startAnnotation) {
        if (nearest_annotation$strand == "+"){
          nearest_annotation$description <- "overlapping 3'"
        } else {
          nearest_annotation$description <- "overlapping 5'"
        }
      } else if (nearest_annotation$end>=nearest_annotation$startAnnotation & 
                 nearest_annotation$end<nearest_annotation$endAnnotation) {
        if (nearest_annotation$strand == "+"){
          nearest_annotation$description <- "overlapping 5'"
        } else {
          nearest_annotation$description <- "overlapping 3'"
        }
      } else {
        nearest_annotation$description <- NA
      }
      nearest_annotation
    }
  ))
}
