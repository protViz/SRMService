context("test-limmahelpers")

test_that("checkSRMService", {
  library(SRMService)
  grp <- SRMService::mqQuantMatrixGRP2
  grp <- mqQuantMatrixGRP2
  tmp <- grp$getResultTable()
  tmp <- grp$getResultTableWithPseudo()

})
