context("test-limmahelpers")

test_that("multiplication works", {
  library(SRMService)
  grp <- SRMService::mqQuantMatrixGRP2
  grp <- mqQuantMatrixGRP2
  tmp <- grp$getResultTable()
  tmp <- grp$getResultTableWithPseudo()

})
