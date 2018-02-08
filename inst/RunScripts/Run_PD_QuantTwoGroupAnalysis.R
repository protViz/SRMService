library(readr)
library(dplyr)
library(quantable)
library(SRMService)

## specify file names
inputFilesFile <- "20161115_02_G3_InputFiles.txt"
proteinsFile <- "20161115_02_G3_Proteins.txt"
##


rm(list=ls())
inputFiles <- read_tsv( inputFilesFile )
inputFileFix <- fix_PD_inputFiles(inputFiles)

annotation <- data.frame(Raw.file = inputFileFix$Raw.File,
                         Condition = gsub("[0-9]","",quantable::split2table(inputFileFix$Raw.File)[,3]),
                         BioReplicate = paste("X",1:nrow(inputFileFix),sep=""),
                         Run = 1:nrow(inputFileFix),
                         IsotopeLabelType = rep("L",nrow(inputFileFix)),
                         stringsAsFactors = F)



proteins <- read_tsv( proteinsFile )
proteinsFIX <- remap_PD_ProteinTable(proteins,inputFiles)




###################################
### Configuration section
fix(annotation)
resultdir="output"
dir.create(resultdir)


Experimentname = ""
nrNas = sum(!is.na(annotation$Condition)) - 1
nrPeptides = 2
reference=unique(annotation$Condition)[1]
qvalueThreshold = 0.01
qfoldchange =2

write.table(annotation, file=file.path(resultdir, "annotationused.txt"))

####### END of user configuration ##



library(SRMService)
grp2 <- Grp2Analysis(annotation, "test pd import", maxNA=nrNas  , nrPeptides=nrPeptides, reference=reference)
grp2$setProteins(proteinsFIX)
grp2$setQValueThresholds(qvalue = qvalueThreshold,qfoldchange = qfoldchange)

results <- grp2$getResultTableWithPseudo()
write_tsv(results,path= file.path(resultdir,"pValues.tsv"))
rmarkdown::render("Grp2Analysis.Rmd", output_format = bookdown::pdf_document2())




