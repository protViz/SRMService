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

#fix(annotation)

Experimentname = ""
nrNas = sum(!is.na(annotation$Condition)) - 1
nrNas = 5
nrPeptides = 2
reference=unique(annotation$Condition)[1]
reference="WT"
qvalueThreshold = 0.05
qfoldchange =1
write.table(annotation, file=file.path(resultdir, "annotationused.txt"))

####### END of user configuration ##

# source("R/Grp2Analysis.R")
grp2 <- Grp2Analysis(annotation, "Experimentname",
                     maxNA=nrNas,
                     nrPeptides=nrPeptides,
                     reference=reference,
                     numberOfProteinClusters = 20
                     )

grp2$getDesignMatrix()

grp2$setMQProteinGroups(protein)
grp2$setQValueThresholds(qvalue = qvalueThreshold,qfoldchange = qfoldchange)
mqQuantMatrixGRP2 <- grp2

head(mqQuantMatrixGRP2$getModPValuesCI())

usethis::use_data(mqQuantMatrixGRP2, overwrite = TRUE)
#readr::write_tsv(grp2$getResultTable(), path=file.path(resultdir,"pValues.csv"))

## REMOVE TO RENDER
# rmarkdown::render("vignettes/Grp2AnalysisHeatmap3.Rmd",bookdown::pdf_document2(), params=list(grp = grp2))
# rmarkdown::render("vignettes/Grp2Analysis.Rmd",bookdown::pdf_document2(), params=list(grp = grp2))

