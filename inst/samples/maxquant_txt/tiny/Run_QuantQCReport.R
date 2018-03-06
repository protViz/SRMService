rm(list=ls())
library(limma)
library(qvalue)
library(SRMService)


protein <- read.table(("proteinGroups.txt"),
                      sep="\t",
                      stringsAsFactors = F,
                      header=T)


rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )


tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp)-5):length(tmp)],ylim=c(0, length(protein$Peptides)),xlab='nr of proteins with at least # peptides')


grp2 <- QCProteinReport( "mein super project", maxNA=3  , nrPeptides=2)
grp2$setMQProteinGroups(protein)

rmarkdown::render("inst/reports/QCReport.Rmd","pdf_document")


