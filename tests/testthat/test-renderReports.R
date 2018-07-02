context("test-renderReports")
tmp_d <- tempdir()

test_that("render TR_SRM_Summary", {
  packagedir <- path.package("SRMService")
  reportFile <- file.path(packagedir, "inst" , "reports" , "TR_SRM_Summary.Rmd")
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
