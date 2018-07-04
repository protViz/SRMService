context("test-renderReports")
library(SRMService)
tmp_d <- tempdir()

packagedir <- path.package("SRMService")

test_that("render TR_SRM_Summary", {
  reportFile <- file.path(packagedir, "inst" , "reports" , "TR_SRM_Summary.Rmd")
  #reportFile <- file.path(packagedir, "reports" , "TR_SRM_Summary.Rmd")

  tmp_f <- file.path(tmp_d,"TR_SRM_Summary.Rmd")
  tmp <- file.copy(reportFile, tmp_f,overwrite = TRUE)
  expect_true(file.exists(reportFile))
  expect_true(file.exists(tmp_f))
  if(expect_true(tmp)){
    skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label.Type", qValue="Detection.Q.Value")
    skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
    data(skylinePRMSampleData)
    sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
    x <- rmarkdown::render(tmp_f, output_format = "html_document",
                      params = list(data = sample_analysis,
                                    configuration = skylineconfig),envir = new.env())

  }
})

test_that("source Run Tidy analysis", {

  reportFile <- file.path(packagedir, "inst" , "RunScripts" , "Run_TidyAnalysis_Skyline_PRM.R")
  source(reportFile)
  })

