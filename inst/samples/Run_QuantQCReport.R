#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following URL
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
#

rm(list=ls())
library(limma)
library(qvalue)
library(SRMService)


protein <- read.table(("inst/samples/proteinGroups/proteinGroups.txt"),
                      sep="\t",
                      stringsAsFactors = F,
                      header=T)


rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )


tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp)-5):length(tmp)],ylim=c(0, length(protein$Peptides)),xlab='nr of proteins with at least # peptides')


grp2 <- QCProteinReport( "p2244_MilenaS_PN_HvsR", maxNA=3  , nrPeptides=2)

grp2$setMQProteinGroups(protein)

rmarkdown::render("inst/reports/QCReport.Rmd","pdf_document")


