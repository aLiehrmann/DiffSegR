## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  data <- newExperiment(
#    sampleInfo   = sample_info,
#    loci         = data.frame(
#      seqid      = c("ChrC", "Chr1"),
#      chromStart = c(71950, 2421089),
#      chromEnd   = c(78500, 2422066),
#      locusID    = c("psbB_petD", "F24B9")
#    ),
#    referenceCondition = "wt",
#    otherCondition     = "pnp1_1",
#    nbThreads          = 4,
#    nbThreadsByLocus   = 2,
#    coverage           = working_directory
#  )

