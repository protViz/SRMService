#
# If you are going to use results produced by the scripts please do cite:
# SRMSerivce R package
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
# FGCZ

#

library(reshape2)
library(limma)
library(SRMService)


# max na per protein
maxNA <- 3
# min number of peptides per protein
nrPeptides <- 2
# Experiment name (will appear as title in the pdf)
experimentName <- "myTesting"




### Do not edit (except you know what you are doing)
protein <- read.table(
  "proteinGroups.txt",
  sep = "\t",
  stringsAsFactors = F,
  header = T
)

tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp) - 5):length(tmp)], ylim = c(0, length(protein$Peptides)), xlab =
          'nr of proteins with at least # peptides')

grp2 <-
  QCProteinReport(experimentName, maxNA = maxNA  , nrPeptides = nrPeptides)
grp2$setMQProteinGroups(protein)


#rmarkdown::render("QCReport.Rmd","pdf_document")
rmarkdown::render("QCReport.Rmd", bookdown::pdf_document2())
