#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following URL
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
#
rm(list=ls())
library(limma)
library(SRMService)

protein <- system.file("samples/proteinGroups/proteinGroupsPullDown.txt",package = "SRMService")
protein <- readr::read_tsv(protein)
colnames(protein) <- make.names(colnames(protein))


tmp <- strsplit(protein$Majority.protein.IDs,split=" ")
tmp2 <- sapply(tmp, function(x){x[1]})
protein$Majority.protein.IDs <- gsub(">","",tmp2)


rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
condition <- quantable::split2table(rawF)
condition <- paste(condition[,4], condition[,5], sep="_")
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

Experimentname = "p2550"
nrNas = 5
nrPeptides = 2

annotation$Condition[grepl("_w$",annotation$Condition)] <- NA
reference = "pegfp_wo" # unique(annotation$Condition)[3]

#write.table(annotation, file="output/annotationused.txt")


####### END of user configuration ##
grp2 <- Grp2Analysis(annotation, "Experimentname", maxNA=nrNas  , nrPeptides=nrPeptides, reference=reference)
grp2$setMQProteinGroups(protein)
grp2$qfoldchange = 2
grp2$setQValueThresholds(qvalue = 0.01)

grp2PullDownExample <- grp2
usethis::use_data(grp2PullDownExample, overwrite = TRUE)
results <- grp2$getResultTable()
#write.table(results, file=file.path(workdir,"pValues.csv"), quote=FALSE, sep = "\t", col.names=NA)

rmarkdown::render("Grp2Analysis.Rmd", bookdown::pdf_document2())
