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
proteinGroupsFile <- "Generic_QuantMatrix.txt"
# read in protein groups file
protein <- readr::read_tsv(proteinGroupsFile)
colnames(protein) <- make.names(colnames(protein))


# important structure for protein matrix
#colnames(protein)
# "Majority.protein.IDs"   "Intensity.dnmt3b_15_s1" "Intensity.dnmt3b_15_s2" "Intensity.dnmt3l_11_s1" "Intensity.dnmt3l_11_s2" "Intensity.IP_lsh_s1"    "Intensity.IP_lsh_s2"    "Intensity.WT"           "Peptides"               "Fasta.headers"   
# fix column names in order to parse all efficientl
# also add some columns that are not by default present
# "Majority.protein.IDs"
# Peptides
# Fasta.headers
originalHeaders <- colnames(protein)
fakeNrPeptides <- rep(3, nrow(protein))
fakeFastaHeader <- rep("No_Description", nrow(protein))
protein <- cbind(protein,fakeNrPeptides, fakeFastaHeader)


# for project p2482 put the sample annotations from metainfo
newHeaders <- c("Majority.protein.IDs", paste("Intensity.",originalHeaders[2:length(originalHeaders)], sep=""), "Peptides","Fasta.headers")
colnames(protein) <- newHeaders




# all raw files in protein groups (select here for proper 2grp)
#rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )

# how to parse condition from the filenames (separator in fn _ )
condition <- quantable::split2table(rawF)[,3]
#
annotation <- data.frame(Raw.file = rawF,
                         Condition = condition,
                         BioReplicate = paste("X",1:length(condition),sep=""),
                         Run = 1:length(condition),
                         IsotopeLabelType = rep("L",length(condition)),
                         stringsAsFactors = F)


###################################
### Configuration section
resultdir <- "p2482_output_TreatedVSUntreated"
dir.create(resultdir)

# calls up data editor
fix(annotation)

# default settings
Experimentname = "p2482 allTNFTreated vs TNFuntreated"
nrNas = sum(!is.na(annotation$Condition)) - 3
nrPeptides = 2
reference=unique(annotation$Condition)[1]
qvalueThreshold = 0.01
qfoldchange =2

write.table(annotation, file=file.path(resultdir, "annotationused.txt"))

####### END of user configuration ##


# important structure for protein matrix
#colnames(protein)
# "Majority.protein.IDs"   "Intensity.dnmt3b_15_s1" "Intensity.dnmt3b_15_s2" "Intensity.dnmt3l_11_s1" "Intensity.dnmt3l_11_s2" "Intensity.IP_lsh_s1"    "Intensity.IP_lsh_s2"    "Intensity.WT"           "Peptides"               "Fasta.headers"   

# Do the analysis
grp2 <- Grp2Analysis(annotation, Experimentname, maxNA=nrNas  , nrPeptides=nrPeptides, reference=reference)
grp2$setMQProteinGroups(protein)
grp2$setQValueThresholds(qvalue = qvalueThreshold,qfoldchange = qfoldchange)

#write out results and render pdf
readr::write_tsv(x = grp2$getResultTable(), path = file.path(resultdir,"pValues.csv"))
rmarkdown::render("Grp2Analysis.Rmd",bookdown::pdf_document2())

#Move 
outPutFileName <- paste(gsub(x = Experimentname, pattern = " ", replacement = "_"),".pdf",sep="")
myMoveCommand <- paste("mv Grp2Analysis.pdf ",resultdir,"/",outPutFileName,sep="")
system(command = myMoveCommand, intern = FALSE)


