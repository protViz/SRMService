rm(list=ls())

protein <- read.csv("d:/googledrive/DataAnalysis/p2084/215579/proteinGroups.txt",sep="\t",stringsAsFactors = F)
rawF <-gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )

condition <- quantable::split2table(rawF)[,3]
annotation <-data.frame(Raw.file = rawF,
                        Condition = condition,
                        BioReplicate = paste("X",1:length(condition),sep=""),
                        Run = 1:length(condition),
                        IsotopeLabelType = rep("L",length(condition)), stringsAsFactors = F)



source("R/Grp2Analysis.R")
grp2 <- Grp2Analysis(annotation, "p2084", maxNA=8 , nrPeptides=2)
grp2$setMQProteinGroups(protein)
table(grp2$getNrNAs())

rmarkdown::render("inst/reports/Grp2Analysis.Rmd",output_format = "pdf_document", output_file = "Grp2Analysis.pdf")
