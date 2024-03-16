#' Annotate Each Provided Differentialy Expressed Region (DER) with Its
#' Closest Annotation(s)
#'
#' @description 
#' `annotateNearest()` relies on user-specified annotations in GFF or GTF 
#' format to annotate the DERs found during the differential expression analysis. 
#' Seven classes of labels translate the relative positions of the stranded DER 
#' and its closest annotation(s): 
#' \itemize{
#'    \item antisense ; 
#'    \item upstream ; 
#'    \item downstream ; 
#'    \item inside ; 
#'    \item overlapping 3’ ; 
#'    \item overlapping 5’ ; 
#'    \item overlapping 5’, 3’.
#' }
#' For the non-stranded DERs, the two classes of labels are: overlapping and 
#' non-overlapping.
#' 
#' @param features A `GenomicRanges` object augmented with the DER column, which 
#' indicates if the region is differentially expressed, as returned by 
#' DiffSegR::dea().
#' @param annotations A `String`. Path to the annotation file (in GFF or GTF 
#' format).
#' @param outputDirectory A `String`. Path to the output directory.
#' @param annotationsType A vector of `String`. Keep only matching annotations.
#' @param orderBy A `String`. The name of the column on which returned 
#' annotations are sorted.
#' @param select A vector of `String`. Keep only columns matching these 
#' names in the returned `Data.frame`. By default all columns from 
#' mcols(dds) and mcols(annotations) are kept.
#' @param verbose A `Boolean`. Should all the operations performed be displayed ?
#' @returns Annotated differentialy expressed regions as `Data.frame` object.
#' 
#' @export
annotateNearest <- function(
  features, 
  annotations,
  outputDirectory = NULL,
  annotationsType = NULL,
  orderBy         = "start",
  select          = NULL,
  verbose         = TRUE) {
  
  features <- features[features$DER==TRUE,]
  if (length(features)<1) {
    stop("No differential expressed region to annotate.")
  }
  
  ## loading annotations
  if (verbose) {
    message(
      "\n > Loading annotations from ",
      annotations,
      " ..."
    )
  }
  annotations <- rtracklayer::import(annotations)
  if (!is.null(annotationsType)) {
    annotations <- annotations[annotations$type %in% annotationsType,]
  }

  stranded <- as.character(BiocGenerics::strand(features)[1]) %in% c("+","-")

  if (verbose) message("\n > Annotate differential expressed regions ...")

  if (stranded) { 
    
    ## (i) DERs and annotations on same strand
    annotated_DERs <- findNearestAnnotation(
      features    = features,
      annotations = annotations,
      antisense   = FALSE
    )

    ## (ii) overlapping DERs and annotations on opposite strand (antisense)
    annotated_DERs_as <- findNearestAnnotation(
      features    = features,
      annotations = annotations,
      antisense   = TRUE
    )

    ## merge (i) and (ii)
    all_annotated_DERs <- rbind(annotated_DERs, annotated_DERs_as)

    ## labeling the annotations
    all_annotated_DERs <- labeling(
      all_annotated_DERs = all_annotated_DERs
    )
  
  } else {
    
    all_annotated_DERs <- findNearestAnnotation(
      features    = features,
      annotations = annotations,
      stranded    = FALSE
    )

    ## labeling the annotations for non-stranded features
    all_annotated_DERs <- labelingNS(
      all_annotated_DERs = all_annotated_DERs
    )
  }

  ## cast list elements in a string
  all_annotated_DERs <- data.frame(lapply(all_annotated_DERs , function(x) {
    if (methods::is(x,"list")) {
      sapply(x, function(elements) {
        ifelse(length(elements)>0, paste(elements, sep=","), '') 
      })
    } else {
      unlist(x)
    }
  }))

  ## keep only selected columns
  select <- select[select%in%names(all_annotated_DERs)]
  if (!is.null(select)) {
    all_annotated_DERs <- all_annotated_DERs[,select]
  }
  
  ## sort the annotations
  all_annotated_DERs <- all_annotated_DERs[order(all_annotated_DERs[[orderBy]]),]

  ## save the annotations
  if (!is.null(outputDirectory)) {
    if (verbose) {
      message(
        "\n > A report can be found in ", 
        file.path(outputDirectory, "annotated_DERs.txt")
      )
    }
    utils::write.table(
      x         = all_annotated_DERs, 
      file      = file.path(outputDirectory, "annotated_DERs.txt"),
      row.names = FALSE
    )
  }
  all_annotated_DERs
}

findNearestAnnotation <- function(
  features,
  annotations,
  antisense = FALSE,
  stranded  = TRUE) {

  ignore_strand <- ifelse(antisense | !stranded, TRUE, FALSE) 

  ## find nearest neighbor annotations 
  nearest_neighbor_tab <- GenomicRanges::distanceToNearest(
    x             = features, 
    subject       = annotations,
    select        = "all",
    ignore.strand = ignore_strand
  )
  
  ## prevents name collisions between annotations and DERs
  nearest_annotations <- as.data.frame(annotations[nearest_neighbor_tab@to,])
  nearest_annotations$startAnnotation  <- nearest_annotations$start
  nearest_annotations$endAnnotation    <- nearest_annotations$end
  nearest_annotations$strandAnnotation <- nearest_annotations$strand
  nearest_annotations$widthAnnotation  <- nearest_annotations$width
  nearest_annotations <- nearest_annotations[
    ,!names(nearest_annotations)%in%c("start","end","strand","width")
  ]
  
  ## merge annotations and DERs
  annotated_DERs <- cbind(
    as.data.frame(
      features[nearest_neighbor_tab@from,],
      row.names = NULL
    ),
    nearest_annotations,
    distance = GenomicRanges::mcols(nearest_neighbor_tab)$distance
  )
  
  ## if antisense keep only overlapping annotations on the oposite strand 
  if (antisense) {
    annotated_DERs <- annotated_DERs[
      annotated_DERs$strand != annotated_DERs$strandAnnotation & 
        annotated_DERs$distance == 0,
    ]
  }
  annotated_DERs
}

labeling <- function(all_annotated_DERs) {
  ## TODO : improve here 
  do.call(rbind, lapply(1:nrow(all_annotated_DERs), function(i) {
    nearest_annotation <- all_annotated_DERs[i,]
    
    if (nearest_annotation$strand != nearest_annotation$strandAnnotation & 
      nearest_annotation$distance == 0) {
      nearest_annotation$description <- "antisense"
    } else if (nearest_annotation$end<nearest_annotation$startAnnotation) {
      if (nearest_annotation$strand == "+") {
        nearest_annotation$description <- "upstream"
      } else {
        nearest_annotation$description <- "downstream"
      }
    } else if (nearest_annotation$start>nearest_annotation$endAnnotation) {
      if (nearest_annotation$strand == "+") {
        nearest_annotation$description <- "downstream"
      } else {  
        nearest_annotation$description <- "upstream"
      }
    } else if (nearest_annotation$start>=nearest_annotation$startAnnotation & 
      nearest_annotation$end<=nearest_annotation$endAnnotation) {
      nearest_annotation$description <- "inside"
    } else if (nearest_annotation$start<=nearest_annotation$startAnnotation & 
      nearest_annotation$end>=nearest_annotation$endAnnotation) {
      nearest_annotation$description <- "overlapping 5', 3'"
    } else if (nearest_annotation$start<=nearest_annotation$endAnnotation & 
      nearest_annotation$start>nearest_annotation$startAnnotation) {
      if (nearest_annotation$strand == "+") {
        nearest_annotation$description <- "overlapping 3'"
      } else {
        nearest_annotation$description <- "overlapping 5'"
      }
    } else if (nearest_annotation$end>=nearest_annotation$startAnnotation & 
      nearest_annotation$end<nearest_annotation$endAnnotation) {
      if (nearest_annotation$strand == "+") {
        nearest_annotation$description <- "overlapping 5'"
      } else {
        nearest_annotation$description <- "overlapping 3'"
      }
    } else {
      nearest_annotation$description <- NA
    }
    nearest_annotation
  }))
}


labelingNS <- function(all_annotated_DERs) {
  ## TODO : improve here 
  do.call(rbind, lapply(1:nrow(all_annotated_DERs), function(i) {
    nearest_annotation <- all_annotated_DERs[i,]
    if (nearest_annotation$distance == 0) {
      nearest_annotation$description <- "overlapping"
    } else {
      nearest_annotation$description <- "non-overlapping"
    }
    nearest_annotation
  }))
}