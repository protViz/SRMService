rm(list=ls())
library(limma)
library(qvalue)
dir("inst/samples/")

#source("R/Grp2Analysis.R")

protein <- read.table(("inst/samples/proteinGroups/proteinGroups2x4.txt"),
                      sep="\t",
                      stringsAsFactors = F,
                      header=T)


grep("Intensity\\.",colnames(protein),value=T)
rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
head(rawF)

source("R/QCReport.R")
head(annotation)

tmp <- cumsum(rev(table(protein$Peptides)))

barplot(tmp[(length(tmp)-5):length(tmp)],ylim=c(0, length(protein$Peptides)),xlab='nr of proteins with at least # peptides')

grp2 <- QCProteinReport( "p2244_MilenaS_PN_HvsR", maxNA=3  , nrPeptides=2)



grp2$setMQProteinGroups(protein)

rmarkdown::render("inst/reports/QCReport.Rmd",bookdown::pdf_document2())

