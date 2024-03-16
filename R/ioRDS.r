exportAsRDS <- function(
  pathToCoverages,
  coverageBySample) {
  for (i_coverage in seq_along(pathToCoverages)) {
    saveRDS(coverageBySample[[i_coverage]], pathToCoverages[[i_coverage]])
  }
}

loadRDS <- function(pathToCoverages) {
  
  coverage_by_sample <- lapply(
    seq_along(pathToCoverages),
    function(i_coverage_path) {
      readRDS(pathToCoverages[[i_coverage_path]])
    }
  )

  coverage_by_sample
}

formatCoverageList <- function(
  sampleInfo,
  coverageBySample) {

  samples                   <- sampleInfo$sample
  names(coverageBySample)   <- samples
  strands                   <- names(coverageBySample[[1]])
  coverage_by_strand        <- vector("list", length(strands))
  names(coverage_by_strand) <- strands 
  
  for (strand in strands) { 
    coverage_by_strand[[strand]]        <- vector("list", length(samples))
    names(coverage_by_strand[[strand]]) <- samples
    for (sample in samples) {
      coverage_by_strand[[strand]][[sample]] <- 
        coverageBySample[[sample]][[strand]]
    }
  }

  coverage_by_strand
}
