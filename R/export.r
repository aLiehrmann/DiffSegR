#' Export DiffSegR results
#'
#' @description 
#' `export()` saves for both strands DERs, not-DERs, segmentation, 
#' coverage profiles and log2-FC per-base in formats readable (bedGraph,gff3) 
#' by genome viewers like the Integrative Genome Viewer (IGV). For IGV, export 
#' also creates an IGV sessiong in xml format that allows loading all tracks in
#' one click.
#' 
#' @param data The `List` object returned by [DiffSegR::loadData()]
#' @param dds The `DESeqDataSet` object returned by [DiffSegR::dea()]. 
#' @param outputDirectory A `String`. Path to the output directory.
#' @param genome A `String`. Path to the genome on which reads have been aligned.
#' @param annotations A `String`. Path to the annotation file.
#' @param igvVersion A `String`. What IGV version should be used in the xml 
#' file. 
#' @param select A vector of `String`. Additional information that appear for 
#' each region in IGV. only columns from mcols(dds) matching these names are
#' kept.
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
#'   sampleInfo   = file.path(working_directory,"sample_info.txt"),
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
#' # Visualizing the DERs.
#' exons_path <- system.file(
#'   "extdata", 
#'   "TAIR10_ChrC_exon.gff3", 
#'   package = "DiffSegR"
#' )
#' ChrC_path  <- system.file(
#'   "extdata", 
#'   "TAIR10_ChrC.fa", 
#'   package = "DiffSegR"
#')
#' export(
#'   data            = data,
#'   dds             = dds, 
#'   outputDirectory = working_directory,
#'   genome          = ChrC_path,
#'   annotations     = exons_path,
#'   select          = c("log2FoldChange", "pvalue", "padj")
#' ) 
#' 
#' # delete working directory 
#' unlink(working_directory, recursive = TRUE)
#' 
#' @export
export <- function(
  data, 
  dds,
  outputDirectory, 
  genome,
  annotations = NA,
  igvVersion = "8",
  select     = c()){
  
  library(rtracklayer)
  outputDirectory  <- normalizePath(outputDirectory)
  select           <- c(select, 'type', 'color', 'ID')
  
  GenomicRanges::mcols(dds)$type <- ifelse(
    GenomicRanges::mcols(dds)$rejectedHypotheses, 
    "DER", 
    "not-DER"
  )
  GenomicRanges::mcols(dds)$color <- ifelse(
    GenomicRanges::mcols(dds)$rejectedHypotheses, 
    ifelse(
      GenomicRanges::mcols(dds)$log2FoldChange>0,
      "85,166,81",
      "154,78,191"
    ),
    "150,150,150"
  )
  GenomicRanges::mcols(dds)$ID <- paste( 
    GenomicRanges::start(dds), 
    GenomicRanges::end(dds),
    sep="_"
  )
  reference      <- data$referenceCondition
  d              <- list(plus="+", minus="-", all="*")
  size_factors   <- DESeq2::sizeFactors(dds)
  i_reference    <- data$sampleInfo$condition == reference
  samples        <- data$sampleInfo$sample
  conditions     <- data$sampleInfo$condition
  chrom          <- as.character(GenomeInfoDb::seqnames(data$locus))
  min_log2_cov_l <- list()
  max_log2_cov_l <- list()
  min_log2_fc_l  <- list()
  max_log2_fc_l  <- list()

  for (strand in names(data$coverages)) { 
    ##- pre-compute scaled coverages by sample -------------------------------## 
    cov <- list()
    for (i in seq_along(samples)) {
      cov[[samples[[i]]]] <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges   = IRanges::IRanges(
	        end    = cumsum(data$coverages[[strand]][[samples[[i]]]]@lengths) +
	                 (GenomicRanges::start(data$locus)-1),
	        width  = data$coverages[[strand]][[samples[[i]]]]@lengths
        ),
        score    = data$coverages[[strand]][[samples[[i]]]]@values/size_factors[[i]]
      )
    }  
    ##- looking for common range to set tracks in igv ------------------------##
    min_cov <- 0
    max_cov <- max(sapply(cov, function(x) max(x$score)))
    min_log2_cov <- 0
    max_log2_cov <- log2(max(sapply(cov, function(x) max(x$score)))+1)
    ##- save coverages and transformed coverages as bedgraph file ------------## 
    for (i in seq_along(samples)) {
      if (conditions[[i]] == reference) {
      	color <- as.integer(c(0,0,255))
      } else { 
        color <- as.integer(c(255,0,0))
      } 
      myTrackLine <- methods::new(
        "GraphTrackLine", 
         type       = "bedGraph", 
         name       = paste(samples[[i]], strand),
         graphType  = "line",
         viewLimits = c(min_cov,max_cov),
         color      = color,
         altColor   = color
      )
      rtracklayer::export.bedGraph(
        cov[[samples[[i]]]], 
        file.path(
          outputDirectory, 
          paste0(samples[[i]],"_", strand,".bedGraph")
        ), 
        trackLine=myTrackLine
      )
      cov[[samples[[i]]]]$score <- log2(cov[[samples[[i]]]]$score+1)
      myTrackLine <- methods::new(
        "GraphTrackLine", 
        type="bedGraph", 
         name=paste0("log2(",samples[[i]],") ",strand),
         graphType  = "line",
         viewLimits = c(min_log2_cov,max_log2_cov),
         color      = color,
         altColor   = color
      )
      rtracklayer::export.bedGraph(
        cov[[samples[[i]]]], 
        file.path(
          outputDirectory, 
          paste0("log2_",samples[[i]],"_",strand,".bedGraph")
        ),
        trackLine=myTrackLine
      )
    }
    ##- mean scaled coverage per condition -----------------------------------##   
    cov <- do.call(cbind, lapply(data$coverages[[strand]], FUN=as.vector))
    scaled_cov <- cov%*%diag(1/size_factors)
    
    avg_ref_rle <- S4Vectors::Rle(
      rowMeans(scaled_cov[,i_reference])
    )
    avg_ref <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(
	      end    = cumsum(avg_ref_rle@lengths) +
	               (GenomicRanges::start(data$locus)-1),
	      width  = avg_ref_rle@lengths
      ),
      score    = avg_ref_rle@values
    )
    other     <- data$sampleInfo$condition[!i_reference][[1]]
    avg_other_rle <- S4Vectors::Rle(
      rowMeans(scaled_cov[,!i_reference])
    )
    avg_other <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(
	      end    = cumsum(avg_other_rle@lengths) + 
	               (GenomicRanges::start(data$locus)-1),
	      width  = avg_other_rle@lengths
      ),
      score    = avg_other_rle@values
    )
    min_cov <- 0
    max_cov <- max(c(avg_other$score, avg_ref$score))
    ##- save as bedGraph file ------------------------------------------------##
    myTrackLine <- methods::new(
      "GraphTrackLine", 
      type       = "bedGraph", 
      name       = paste(reference,strand),
      graphType  = "line",
      viewLimits = c(min_cov,max_cov),
      color      = as.integer(c(0,0,255)),
      altColor   = as.integer(c(0,0,255))
    )
    rtracklayer::export.bedGraph(
      avg_ref, 
      file.path(
        outputDirectory, 
        paste0(reference,"_",strand,".bedGraph")
      ), 
      trackLine=myTrackLine
    )
    myTrackLine <- methods::new(
      "GraphTrackLine", 
      type="bedGraph", 
      name=paste(other,strand),
      graphType  = "line",
      viewLimits = c(min_cov,max_cov),
      color      = as.integer(c(255,0,0)),
      altColor   = as.integer(c(255,0,0))
    )
    rtracklayer::export.bedGraph(
      avg_other, 
      file.path(
        outputDirectory, 
        paste0(other,"_",strand,".bedGraph")
      ), 
      trackLine=myTrackLine
    )  
    ##- log2 mean scaled coverage per condition ------------------------------##   
    cov <- do.call(cbind, lapply(data$coverages[[strand]], FUN=as.vector))
    scaled_cov <- cov%*%diag(1/size_factors)
    
    avg_ref_rle <- S4Vectors::Rle(
      rowMeans(log2(scaled_cov[,i_reference]+1))
    )
    avg_ref <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(
	      end    = cumsum(avg_ref_rle@lengths) +
	               (GenomicRanges::start(data$locus)-1),
	      width  = avg_ref_rle@lengths
      ),
      score    = avg_ref_rle@values
    )
    other     <- data$sampleInfo$condition[!i_reference][[1]]
    avg_other_rle <- S4Vectors::Rle(
      rowMeans(log2(scaled_cov[,!i_reference]+1))
    )
    avg_other <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(
	      end    = cumsum(avg_other_rle@lengths) +
	               (GenomicRanges::start(data$locus)-1),
	      width  = avg_other_rle@lengths
      ),
      score    = avg_other_rle@values
    )
    min_log2_cov <- 0
    max_log2_cov <- max(c(avg_other$score, avg_ref$score))
    min_log2_cov_l[[strand]] <- as.character(round(min_log2_cov,6)) 
    max_log2_cov_l[[strand]] <- as.character(round(max_log2_cov,6))
    ##- save as bedGraph file ------------------------------------------------##
    myTrackLine <- methods::new(
      "GraphTrackLine", 
      type       = "bedGraph", 
      name       = paste("log2(",reference,")",strand),
      graphType  = "line",
      viewLimits = c(min_log2_cov,max_log2_cov),
      color      = as.integer(c(0,0,255)),
      altColor   = as.integer(c(0,0,255))
    )
    rtracklayer::export.bedGraph(
      avg_ref, 
      file.path(
        outputDirectory, 
        paste0("log2_",reference,"_",strand,".bedGraph")
      ), 
      trackLine=myTrackLine
    )
    myTrackLine <- methods::new(
      "GraphTrackLine", 
      type="bedGraph", 
      name=paste("log2(",other,")",strand),
      graphType  = "line",
      viewLimits = c(min_log2_cov,max_log2_cov),
      color      = as.integer(c(255,0,0)),
      altColor   = as.integer(c(255,0,0))
    )
    rtracklayer::export.bedGraph(
      avg_other, 
      file.path(
        outputDirectory, 
        paste0("log2_",other,"_",strand,".bedGraph")
      ), 
      trackLine=myTrackLine
    )   
    ##- log2-FC --------------------------------------------------------------##
    scaled_cov2 <- cov%*%diag(1/size_factors)
    log2_fc_rle <- S4Vectors::Rle(
      rowMeans(log2(scaled_cov2[,!i_reference]+1)) - 
      rowMeans(log2(scaled_cov2[,i_reference]+1))
    )
    log2_fc <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(
	      end    = cumsum(log2_fc_rle@lengths) +
	               (GenomicRanges::start(data$locus)-1),
	      width  = log2_fc_rle@lengths
      ),
      score    = log2_fc_rle@values
    )
    ##- save as bedGraph file ------------------------------------------------## 
    myTrackLine <- methods::new(
      "GraphTrackLine", 
      type="bedGraph", 
      name=paste("log2-FC",strand),
      graphType  = "line",
      color      = as.integer(c(149,149,149)),
      altColor   = as.integer(c(149,149,149))
    )
    rtracklayer::export.bedGraph(
      log2_fc, 
      file.path(
        outputDirectory, 
        paste0("log2_FC_",strand,".bedGraph")
      ), 
      trackLine=myTrackLine
    )
    ##- denoised log2-FC -----------------------------------------------------##
    model <- SummarizedExperiment::rowRanges(
      dds[GenomicRanges::strand(dds)==d[[strand]],]
    )
    GenomicRanges::mcols(model)$score <- 
    	GenomicRanges::mcols(model)$log2FoldChangeMean
    ##- save as bedGraph file ------------------------------------------------##
    myTrackLine <- methods::new(
      "GraphTrackLine", 
      type="bedGraph", 
      name=paste("model",strand),
      graphType  = "line",
      viewLimits = c(min(log2_fc$score),max(log2_fc$score)),
      color      = as.integer(c(0,21,150)),
      altColor   = as.integer(c(0,21,150))
    )
    rtracklayer::export.bedGraph(
      model, 
      file.path(
        outputDirectory, 
        paste0("model_",strand,".bedGraph")
      ), 
      trackLine=myTrackLine
    )
    min_log2_fc_l[[strand]] <- as.character(round(min(log2_fc$score),6))
    max_log2_fc_l[[strand]] <- as.character(round(max(log2_fc$score),6))
    ##- save features as gff3 file -------------------------------------------##
    if (length(select)==3) {
      rtracklayer::export(
        SummarizedExperiment::rowRanges(dds[GenomicRanges::strand(dds)==d[[strand]],]), 
        file.path(
          outputDirectory, 
          paste0("DiffSegR_results_", strand,".gff3")
        ),
        "gff3"
      )
    } else {
      rtracklayer::export(
        SummarizedExperiment::rowRanges(dds[GenomicRanges::strand(dds)==d[[strand]],])[,
        select], 
        file.path(
          outputDirectory, 
          paste0("DiffSegR_results_", strand,".gff3")
        ),
        "gff3"
      )
    }
  }
  ##- create igv session from blueprint --------------------------------------##
  igv_session_bp <- xml2::read_xml(
    system.file("extdata", "igv_session_blueprint.xml", package = "DiffSegR")
  )
  xml2::xml_set_attr(igv_session_bp, attr="genome", value=genome)
  xml2::xml_set_attr(igv_session_bp, attr="locus", value= paste0(
      as.character(GenomeInfoDb::seqnames(data$locus)),
      ":",
      GenomicRanges::start(data$locus),
      "-",
      GenomicRanges::end(data$locus)
    )
  )
  xml2::xml_set_attr(igv_session_bp, attr="path", value=file.path(
    outputDirectory,
    "DiffSegR_results_igv_session.xml"
  ))
  xml2::xml_set_attr(igv_session_bp, attr="version", value=igvVersion)
  xml2::xml_set_attr(
    xml2::xml_children(xml2::xml_children(igv_session_bp)[1]), 
    attr = "path",
    value = file.path(
      outputDirectory,
      c(
        paste0("log2_",other,"_plus.bedGraph"),
        "DiffSegR_results_minus.gff3",
        "model_minus.bedGraph",
        "log2_FC_minus.bedGraph",
        paste0("log2_",reference,"_minus.bedGraph"),
        paste0("log2_",reference,"_plus.bedGraph"),    
        paste0("log2_",other,"_minus.bedGraph"),
        "log2_FC_plus.bedGraph",      
        "DiffSegR_results_plus.gff3",
        annotations,     
        "model_plus.bedGraph"
      )
    )
  )
  xml2::xml_set_attr(
    xml2::xml_children(xml2::xml_children(igv_session_bp)[1])[10], 
    attr = "path",
    value = annotations
  )
  xml2::xml_set_attr(
    xml2::xml_children(xml2::xml_children(xml2::xml_children(igv_session_bp))), 
    attr="id",
    value = file.path(
      outputDirectory,
      c(
         paste0("log2_",reference,"_plus.bedGraph"),
         paste0("log2_",other,"_plus.bedGraph"),
         "log2_FC_plus.bedGraph",
         "model_plus.bedGraph",
         "model_minus.bedGraph",
         "log2_FC_minus.bedGraph",
         paste0("log2_",other,"_minus.bedGraph"),      
         paste0("log2_",reference,"_minus.bedGraph")
       )
     )
   )
  
  xml2::xml_set_attr(
    xml2::xml_children(xml2::xml_children(igv_session_bp))[c(14,16)], 
    att="id",
    value = file.path(
      outputDirectory,
      xml2::xml_attr(xml2::xml_children(xml2::xml_children(igv_session_bp))[c(14,16)], att="id") 
    )
  )
  xml2::xml_set_attr(
    xml2::xml_children(xml2::xml_children(igv_session_bp))[15], 
    att="id",
    value = annotations
  )
  xml2::write_xml(
    igv_session_bp, 
    file.path(outputDirectory, "DiffSegR_results_igv_session.xml")
  )
}
