rm(list=ls())
library(limma)
library(qvalue)
library(SRMService)

protein <- read.table(("inst/samples/proteinGroups/proteinGroups2x4.txt"),
                     sep="\t",
                     stringsAsFactors = F,
                     header=T)

# protein <- read.table(("d:/projects/p2244_MilenaS_PN/Claudia20170313/test all_H_R.txt"),
#                       sep="\t",
#                       stringsAsFactors = F,
#                       header=T)

rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
condition <- quantable::split2table(rawF)[,3]
annotation <-data.frame(Raw.file = rawF,
                        Condition = condition,
                        BioReplicate = paste("X",1:length(condition),sep=""),
                        Run = 1:length(condition),
                        IsotopeLabelType = rep("L",length(condition)),
                        stringsAsFactors = F)

# If you want to change the labelling
# fix(annotation)

tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp)-5):length(tmp)],ylim=c(0, length(protein$Peptides)),xlab='nr of proteins with at least # peptides')

grp2 <- Grp2Analysis(annotation, "p2244_MilenaS_PN_HvsR", maxNA=3  , nrPeptides=2, reference="WT")
grp2$setMQProteinGroups(protein)

rmarkdown::render("inst/reports/Grp2Analysis.Rmd",bookdown::pdf_document2())

