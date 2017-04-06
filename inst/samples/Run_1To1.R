##
##

maxquanttxtdirectory <- ""
evidence <- read.table("evidence.txt")

rmarkdown::render("inst/reports/QC1To1.Rmd","pdf_document")
