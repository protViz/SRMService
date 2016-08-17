rm(list=ls())

protein <- read.csv("d:/googledrive/DataAnalysis/p2084/215579/proteinGroups.txt",sep="\t",stringsAsFactors = F)

colnames(protein)

rawF <-gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )

condition <- quantable::split2table(rawF)[,3]
annotation <-data.frame(Raw.file = rawF,
                        Condition = condition,
                        BioReplicate = paste("X",1:length(condition),sep=""),
                        Run = 1:length(condition),
                        IsotopeLabelType = rep("L",length(condition)), stringsAsFactors = F)


library(SRMService)

grp2 <- Grp2Analysis(annotation, "p2084BlaBla", maxNA=8  , nrPeptides=2)
grp2$setMQProteinGroups(protein)

library(dScipa)
rmarkdown::render("inst/reports/Grp2Analysis.Rmd",output_format = "pdf_document", output_file = "Grp2Analysis.pdf")

