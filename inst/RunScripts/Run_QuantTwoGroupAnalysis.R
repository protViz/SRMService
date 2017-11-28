#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following URL
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
#
rm(list=ls())
library(limma)
library(SRMService)

protein <- read.table(("proteinGroups.txt"),
                      sep="\t",
                      stringsAsFactors = F,
                      header=T)

rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
condition <- quantable::split2table(rawF)[,3]
annotation <- data.frame(Raw.file = rawF,
                         Condition = condition,
                         BioReplicate = paste("X",1:length(condition),sep=""),
                         Run = 1:length(condition),
                         IsotopeLabelType = rep("L",length(condition)),
                         stringsAsFactors = F)

workdir <- "output"
dir.create(workdir)

tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp)-5):length(tmp)],ylim=c(0, length(protein$Peptides)),xlab='nr of proteins with at least # peptides')

###################################
### Configuration section

Experimentname = ""
nrNas = 3
nrPeptides = 2
reference=unique(annotation$Condition)[1]
fix(annotation)

write.table(annotation, file="output/annotationused.txt")

####### END of user configuration ##


grp2 <- Grp2Analysis(annotation, "Experimentname", maxNA=nrNas  , nrPeptides=nrPeptides, reference=reference)
grp2$setMQProteinGroups(protein)

write.table(grp2$getResultTable(), file=file.path(workdir,"pValues.csv"), quote=FALSE, sep = "\t", col.names=NA)

rmarkdown::render("Grp2Analysis.Rmd",bookdown::pdf_document2())

