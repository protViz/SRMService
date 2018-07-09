context("test-renderdiareports")

packagedir <- path.package("SRMService")
test_that("source Run DIA correlation analysis", {
  reportFile <- file.path(packagedir, "inst" , "RunScripts" , "Run_TidyCorrelationAnalysis.R")
  source(reportFile)
})
