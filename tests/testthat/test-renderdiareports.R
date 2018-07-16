context("test-renderdiareports")

packagedir <- path.package("SRMService")
#packagedir <- system.file(package = "SRMService")
#packagedir = "."


test_that("run MSqRob", {
  if(require("MSqRob")){
  reportFile <- file.path(packagedir, "inst" , "RunScripts" , "Run_Tidy_MSqRob_Analysis.R")
  source(reportFile)
  }
})
