#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following URL
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
#
rm(list=ls())
library(limma)
library(SRMService)

### Protein groups file
packagedir <- path.package("SRMService")

proteinGroupsFile <- file.path(packagedir, "samples/proteinGroups/proteinGroups.txt")

###


protein <- readr::read_tsv(proteinGroupsFile)
colnames(protein) <- make.names(colnames(protein))
tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp)-5):length(tmp)],ylim=c(0, length(protein$Peptides)),xlab='nr of proteins with at least # peptides')


rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
condition <- quantable::split2table(rawF)[,3]
annotation <- data.frame(Raw.file = rawF,
                         Condition = condition,
                         BioReplicate = paste("X",1:length(condition),sep=""),
                         Run = 1:length(condition),
                         IsotopeLabelType = rep("L",length(condition)),
                         stringsAsFactors = F)



###################################
### Configuration section
resultdir <- "output"
dir.create(resultdir)

fix(annotation)


Experimentname = ""
nrNas = sum(!is.na(annotation$Condition)) - 1
nrNas = 5
nrPeptides = 2
reference=unique(annotation$Condition)[1]
reference="NT"
qvalueThreshold = 0.05
qfoldchange =1

write.table(annotation, file=file.path(resultdir, "annotationused.txt"))

####### END of user configuration ##


grp2 <- Grp2Analysis(annotation, "Experimentname", maxNA=nrNas  , nrPeptides=nrPeptides, reference=reference)
grp2$setMQProteinGroups(protein)
grp2$setQValueThresholds(qvalue = qvalueThreshold,qfoldchange = qfoldchange)

readr::write_tsv(grp2$getResultTable(), path=file.path(resultdir,"pValues.csv"))

rmarkdown::render("Grp2Analysis.Rmd",bookdown::pdf_document2())

