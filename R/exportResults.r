#' Export DiffSegR Results
#'
#' @description 
#' For each locus, `exportResults()` saves DERs, not-DERs, segmentation, 
#' coverage profiles and Log2-FC per-base in formats readable (bedGraph,GFF3) 
#' by genome viewers like the Integrative Genome Viewer (IGV). For IGV, export 
#' also creates an IGV sessiong in XML that allows loading all tracks in
#' one click. exportResults() can operate on each locus in parallel.
#' 
#' @param data The `List` object returned by [DiffSegR::newExperiment()].
#' @param features The `GRanges` part of the object returned by 
#' [DiffSegR::dea()]. 
#' @param outputDirectory A `String`. Path to the output directory. One directory 
#' by locus will be create within outputDirectory.
#' @param genome A `String`. Path to the genome on which reads have 
#' been aligned.
#' @param sizeFactors A `numeric`vector of normalization factors for samples.` 
#' @param annotations A `String`. Path to the annotation file.
#' @param igvVersion A `String`. What IGV version should be used in the XML 
#' file. 
#' @param select A vector of `String`. Additional information that appear for 
#' each region in IGV. only columns from features matching these names are
#' kept.
#' @param verbose A `Logical`. Should all the operations performed be 
#' displayed ?
#'
#' @export
exportResults <- function(
  data, 
  features,
  sizeFactors,
  outputDirectory, 
  genome,
  annotations = NA,
  igvVersion = "8",
  select     = c(),
  verbose    = TRUE) {
  
  require(rtracklayer)
  outputDirectory <- normalizePath(outputDirectory)
  genome          <- normalizePath(genome)
  annotations     <- normalizePath(annotations)
  select          <- c(select, 'type', 'color', 'ID')
  dico_strand     <- list(plus="+", minus="-", all="*")

  if (verbose & data$nbThreads %/% data$nbThreadsByLocus>1) {
    
    ## One message for all since the one inside the loop will not be displayed
    message("\n > Export results from all loci ...")
  }

  start_time <- Sys.time() 
  
  tmp <- do.call(c,customLapply(1:length(data$loci), function(i_locus) {
    
    current_locus <- data$loci[i_locus,] 

    if (verbose) {
      message(
        "\n > Export results from locus ",
        as.character(GenomicRanges::seqnames(current_locus)),
        " ",
        GenomicRanges::start(current_locus),
        "-",
        GenomicRanges::end(current_locus),
        " (ID : ",
        current_locus$locusID,
        ") ..."
      )
    }

    current_output_dir <- file.path(outputDirectory, current_locus$locusID)
    if(!file.exists(current_output_dir)) {
      dir.create(file.path(outputDirectory, current_locus$locusID))
    }

    target_samples <- data$sampleInfo[
      data$sampleInfo$condition %in% c(data$referenceCondition, data$otherCondition),]

    i_reference <- target_samples$condition == data$referenceCondition

    path_to_coverages <- sub(
      ".rds", 
      paste0("_", current_locus$locusID, ".rds"),
      target_samples$coverage
    )
    
    coverage_by_sample <- loadRDS(pathToCoverages = path_to_coverages)

    coverage_by_strand <- formatCoverageList(
      sampleInfo       = target_samples,
      coverageBySample = coverage_by_sample
    )

    strands <- names(coverage_by_strand)

    for (strand in strands) {
    
      ## estimate scaled/normalized coverages 
      scaled_coverages <- lapply(seq_along(coverage_by_strand[[strand]]), function(sample) {
        coverage_by_strand[[strand]][[sample]]/sizeFactors[[sample]]
      })

      ## export scaled mean coverage by condition in bedgraph files 
      exportAvgLog2Coverage(
        locus            = current_locus,
        scaled_coverages = scaled_coverages,
        strand           = strand,
        ref              = data$referenceCondition,
        other            = data$otherCondition,
        i_reference      = i_reference,
        outputDirectory  = current_output_dir
      )

      target_features <- features[GenomicRanges::findOverlaps(current_locus, features)@to,]

      ## export scaled log2-FC and segmentation model in bedgraph files  
      exportLog2FoldChange(
        locus            = current_locus,
        features         = target_features[GenomicRanges::strand(target_features)==dico_strand[[strand]],],
        scaled_coverages = scaled_coverages,
        i_reference      = i_reference,
        strand           = strand,
        dicoStrand       = dico_strand,
        outputDirectory  = current_output_dir
      )

      ## export up-regulated / down-regulated and not-diff regions in bed graph file
      exportDeaResults(
        features        = target_features[GenomicRanges::strand(target_features)==dico_strand[[strand]],],
        select          = select,
        strand          = strand,
        outputDirectory = current_output_dir
      )
    }

    ## create an xml igv session 
    exportIgvSession(
      locus           = current_locus,
      genome          = genome,
      annotations     = annotations,
      ref             = data$referenceCondition,
      other           = data$otherCondition,
      igvVersion      = igvVersion,
      outputDirectory = current_output_dir,
      stranded        = ifelse(strands[[1]]=="all", FALSE, TRUE)
    )
  }, nbThreads = data$nbThreads %/% data$nbThreadsByLocus))
}  


exportAvgLog2Coverage <- function(
  locus,
  scaled_coverages,
  i_reference,
  strand,
  ref,
  other,
  outputDirectory) {
  
  ## estimate the averaged scaled or normalized coverages across replicates
  ## reference condition
  scaled_log2_ref   <- transformationFactory(type = "avg.Rle")(
    scaled_coverages[i_reference]
  )

  ## other condition
  scaled_log2_other <- transformationFactory(type = "avg.Rle")(
    scaled_coverages[!i_reference]
  )

  ## cast as features objects
  ## reference condition
  scaled_log2_ref   <- rleAsfeatures(scaled_log2_ref, locus)

  ## other condition
  scaled_log2_other <- rleAsfeatures(scaled_log2_other, locus)

  ## this will allow for the fixing of a common viewLimits issue in IGV tracks
  min_scaled_log2_cov <- 0
  max_scaled_log2_cov <- max(scaled_log2_other$score, scaled_log2_ref$score)
  
  ## save as bedGraph file
  ## reference condition
  myTrackLine <- methods::new(
    "GraphTrackLine", 
    type       = "bedGraph", 
    name       = paste("log2(",ref,")",strand),
    graphType  = "line",
    viewLimits = c(min_scaled_log2_cov, max_scaled_log2_cov),
    color      = as.integer(c(0,0,255)),
    altColor   = as.integer(c(0,0,255))
  )
  rtracklayer::export.bedGraph(
    scaled_log2_ref, 
    file.path(
      outputDirectory, 
      paste0("log2_",ref,"_",strand,".bedGraph")
    ), 
    trackLine = myTrackLine
  )
  
  ## other condition
  myTrackLine <- methods::new(
    "GraphTrackLine", 
    type       = "bedGraph", 
    name       = paste("log2(" ,other, ")", strand),
    graphType  = "line",
    viewLimits = c(min_scaled_log2_cov, max_scaled_log2_cov),
    color      = as.integer(c(255,0,0)),
    altColor   = as.integer(c(255,0,0))
  )
  rtracklayer::export.bedGraph(
    scaled_log2_other, 
    file.path(
      outputDirectory, 
      paste0("log2_",other,"_",strand,".bedGraph")
    ), 
    trackLine = myTrackLine
  )
}

exportLog2FoldChange <- function(
  locus,
  features,
  scaled_coverages,
  i_reference,
  strand,
  dicoStrand,
  outputDirectory) {

  ## estimate scaled/normalized log2-FC per-base 
  scaled_lfc <- transformationFactory(type = "log2FoldChange.Rle")(
    numerator   = scaled_coverages[!i_reference],
    denominator = scaled_coverages[i_reference]
  )

  ## estimate denoised log2-FC
  starts <- features$modelStart
  ends   <- features$modelEnd
  lfc_mean <- sapply(
    seq_along(starts),
    function(i) mean(scaled_lfc[starts[[i]]:ends[[i]]])
  )
  features$score <- lfc_mean

  ## Find gaps that can appear after removing regions with low expression (for 
  ## instance), and fill them.
  ## NOTE : Gaps are already filled at genomeAnalysis level. 
  locus_tmp <- locus
  BiocGenerics::strand(locus_tmp) <- dicoStrand[[strand]]
  complementary <- GenomicRanges::setdiff(locus_tmp, features)
  if(length(complementary)>0) {
    complementary$score<-0
    features <- sort(c(features, complementary))
  }

  ## save as bedGraph file
  myTrackLine <- methods::new(
    "GraphTrackLine", 
    type       = "bedGraph", 
    name       = paste("model", strand),
    graphType  = "line",
    viewLimits = c(min(scaled_lfc), max(scaled_lfc)),
    color      = as.integer(c(0,21,150)),
    altColor   = as.integer(c(0,21,150))
  )
  rtracklayer::export.bedGraph(
    features, 
    file.path(
      outputDirectory, 
      paste0("model_", strand, ".bedGraph")
    ), 
    trackLine = myTrackLine
  )

  ## save the log2-FC profile as bedgraph file 
  scaled_lfc <- rleAsfeatures(scaled_lfc, locus)
 
  myTrackLine <- methods::new(
    "GraphTrackLine", 
    type       = "bedGraph", 
    name       = paste("log2-FC", strand),
    graphType  = "line",
    viewLimits = c(min(scaled_lfc$score), max(scaled_lfc$score)),
    color      = as.integer(c(149,149,149)),
    altColor   = as.integer(c(149,149,149))
  )
  rtracklayer::export.bedGraph(
    scaled_lfc, 
    file.path(
      outputDirectory, 
      paste0("log2_FC_", strand, ".bedGraph")
    ), 
    trackLine = myTrackLine
  )
}


exportDeaResults <- function(
  features,
  select,
  strand,
  outputDirectory) {

  features$type <- ifelse(
    features$DER, 
    "DER", 
    "not-DER"
  )
  features$color <- ifelse(
    features$DER, 
    ifelse(
      features$log2FoldChange>0,
      "85,166,81", ## up-regulated : green 
      "154,78,191" ## down-regulated : purple
    ),
    "150,150,150" ## not-diff : grey
  )
  features$ID <- paste(
    GenomicRanges::start(features), 
    GenomicRanges::end(features),
    sep="_"
  )

  if (length(select)==3) {  ## by default select = c('type', 'color', 'ID') 
    
    ## export found regions and annotate them with every informations from features 
    rtracklayer::export(
      features, 
      file.path(
        outputDirectory, 
        paste0("DiffSegR_results_", strand,".gff3")
      ),
      "gff3"
    )
  } else {

    ## export found regions and annotate them with only user-selected informations
    ## from features
    rtracklayer::export(
      features[,select], 
      file.path(
        outputDirectory, 
        paste0("DiffSegR_results_", strand,".gff3")
      ),
      "gff3"
    )
  }
}


exportIgvSession <- function(
  locus,
  genome,
  annotations,
  ref,
  other,
  igvVersion,
  outputDirectory,
  stranded) {
  
  if (stranded) {
    
    igv_session_bp <- xml2::read_xml(system.file("extdata", "igv_session_blueprint.xml", package = "DiffSegR"))
    xml2::xml_set_attr(igv_session_bp, attr = "genome", value = genome)
    xml2::xml_set_attr(igv_session_bp, attr = "locus", value = paste0(as.character(GenomeInfoDb::seqnames(locus)), ":", GenomicRanges::start(locus), "-", GenomicRanges::end(locus)))
    xml2::xml_set_attr(igv_session_bp, attr = "path", value = file.path(".", "DiffSegR_results_igv_session.xml"))
    xml2::xml_set_attr(igv_session_bp, attr = "version", value = igvVersion)
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp)[1]), attr = "path", value = file.path(".",
      c(
        paste0("log2_",other,"_plus.bedGraph"),
        "DiffSegR_results_minus.gff3",
        "model_minus.bedGraph",
        "log2_FC_minus.bedGraph",
        paste0("log2_",ref,"_minus.bedGraph"),
        paste0("log2_",ref,"_plus.bedGraph"),    
        paste0("log2_",other,"_minus.bedGraph"),
        "log2_FC_plus.bedGraph",      
        "DiffSegR_results_plus.gff3",
        annotations,     
        "model_plus.bedGraph"
      )
    ))
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp)[1])[10], attr = "path", value = annotations)
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(xml2::xml_children(igv_session_bp))), attr = "id", value = file.path(".",
      c(
        paste0("log2_",ref,"_plus.bedGraph"),
        paste0("log2_",other,"_plus.bedGraph"),
        "log2_FC_plus.bedGraph",
        "model_plus.bedGraph",
        "model_minus.bedGraph",
        "log2_FC_minus.bedGraph",
        paste0("log2_",other,"_minus.bedGraph"),      
        paste0("log2_",ref,"_minus.bedGraph")
      )
    ))
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[c(14,17)], att = "id", value = file.path(".", xml2::xml_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[c(14,17)], att="id")))
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[15], att = "id", value = genome)
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[16], att = "id", value = annotations)
    xml2::write_xml(igv_session_bp, file.path(outputDirectory, "DiffSegR_results_igv_session.xml"))
  
  } else {

    igv_session_bp <- xml2::read_xml(system.file("extdata", "igv_session_blueprint.xml", package = "DiffSegR"))
    xml2::xml_set_attr(igv_session_bp, attr = "genome", value = genome)
    xml2::xml_set_attr(igv_session_bp, attr = "locus", value = paste0(as.character(GenomeInfoDb::seqnames(locus)), ":", GenomicRanges::start(locus), "-", GenomicRanges::end(locus)))
    xml2::xml_set_attr(igv_session_bp, attr = "path", value = file.path(".", "DiffSegR_results_igv_session.xml"))
    xml2::xml_set_attr(igv_session_bp, attr = "version", value = igvVersion)
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp)[1]), attr = "path", value = file.path(".",
      c(
        paste0("log2_",other,"_all.bedGraph"),
        paste0("log2_",ref,"_all.bedGraph"), 
        "log2_FC_all.bedGraph",      
        "DiffSegR_results_all.gff3",
        annotations,     
        "model_all.bedGraph"
      )
    ))
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp)[1])[5], attr = "path", value = annotations)
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(xml2::xml_children(igv_session_bp))), attr = "id", value = file.path(".",
      c(
        paste0("log2_",ref,"_all.bedGraph"),
        paste0("log2_",other,"_all.bedGraph"),
        "log2_FC_all.bedGraph",
        "model_all.bedGraph"
      )
    ))
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[9], att = "id", value = file.path(".", xml2::xml_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[9], att="id")))
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[10], att = "id", value = genome)
    xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[11], att = "id", value = annotations)
    xml2::write_xml(igv_session_bp, file.path(outputDirectory, "DiffSegR_results_igv_session.xml"))
    
  }
}

rleAsfeatures <- function(rleObj, locus) {
    GenomicRanges::GRanges(
    seqnames = as.character(GenomeInfoDb::seqnames(locus)),
    ranges   = IRanges::IRanges(
	    end    = cumsum(rleObj@lengths)+GenomicRanges::start(locus)-1,
	    width  = rleObj@lengths
    ),
    score    = rleObj@values
  )
}