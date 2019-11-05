# R
# 2019-09-23
# 2019-11-05 -> working locally

#
#           This script should be used for the PhosphoService Analysis
#           -> the workflow integrates the PhosphoSTY and the global input protein file (proteinGroups.txt)
#           -> statistics are done on phospho peptide (and global proteome) level, visualized, pdf-reported
#

rm(list=ls())

#################################################
#     Configuration
#################################################

# Variables
LocProbThreshold <- 0.75
nrPeptides = 2
maxNAsForComparison <- 4
min_nonNAinTotal <- maxNAsForComparison
pModThreshold <- 0.2
qvalueThreshold = 0.2
absLg2FCThreshold <- 0.57

# for reporting
fgczProject <<- "pMyProject"
projectName <<- "PhosphoExperiment"

# files, global filters and paramteters and thresholds
### PhosphoSTY file
myPhosphoSTY = system.file("samples/phospho/Phospho_STY_Sites.txt", package = "SRMService")
myAnno_phosphoFile <- system.file("samples/phospho/phospho_TwoGroup_annotationused.txt", package = "SRMService")
###


### Protein groups file
proteinGroupsFile <- system.file("samples/phospho/proteinGroups.txt", package = "SRMService")
myAnno_proteinFile <- system.file("samples/phospho/protein_TwoGroup_annotationused.txt", package = "SRMService")

###
# read in files
# annotations
myAnno <- read.table(myAnno_phosphoFile, header=TRUE)
myAnno_protein <-  read.table(myAnno_proteinFile, header=TRUE)

# matrices
protein <- readr::read_tsv(proteinGroupsFile)
colnames(protein) <- make.names(colnames(protein))




#
myMat_good <- prepareSmallerMatrixFromPhospohSTY(phosphoSTYtableFromMQ = myPhosphoSTY)

# Run Report
PDFfn <- paste(fgczProject, "_PhosphoPeptideOverview_", projectName, ".pdf", sep="")
rmarkdown::render("fgcz_phosphoPeptides_OverviewReport.Rmd",bookdown::pdf_document2())
file.copy(from = "fgcz_phosphoPeptides_OverviewReport.pdf", to = PDFfn, overwrite=TRUE)
file.remove("fgcz_phosphoPeptides_OverviewReport.pdf")






#################################################
#          Quantitative testing - PhosphoPeptideLevel! --> only 1 combination
#################################################
myQuantMatForTesting <- prepareGenericIntensityMatrixFromGoodMatrixForTesting(myMat_good)

# write out matrix to read it back in (hack)
write.table(x = myQuantMatForTesting, file = "pXXXX_FullQuantMatrixUsedForAllTesting.txt", sep="\t", row.names = FALSE,
            quote = FALSE)

# read matrix back in (less problem) -> Error in strsplit(MQProteinGroups$Majority.protein.IDs, split = ";") : at 2grp comparison
phosphopepMatrix <- readr::read_tsv("pXXXX_FullQuantMatrixUsedForAllTesting.txt")
colnames(phosphopepMatrix) <- make.names(colnames(phosphopepMatrix))

# 2 grp comparison
(myRef <- unique(myAnno$Condition)[3])
quantResPhosphoPep <- get2GroupAnalysisFromExistingAnnotationFileNoRender(proteinMat = phosphopepMatrix,
                                                                                            expName = "KOvsCtrl",
                                                                                            resultDirectory = "pXXXX_KOvsCtrl",
                                                                                            pathNFileForAnnotation = myAnno_phosphoFile,
                                                                                            mxNAs = maxNAsForComparison, nrPep = nrPeptides,
                                                                                            DenomReference = as.character(myRef),
                                                                                            lg2FCthreshold = absLg2FCThreshold,
                                                                                            qvThreshold = qvalueThreshold,
                                                                                            Pairing = FALSE)

#################################################
#                 GLOBAL PROTEIN QUANTIFICATION -> Quantitative testing - Global Proteome Level
#################################################

# work on annotations
# rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
#
# myGlobal_condition <- quantable::split2table(rawF)
#
# parsedCondition <- paste(myGlobal_condition[,1],myGlobal_condition[,2], sep="_")
#
# Full_annotation_protein <- data.frame(Raw.file = rawF,
#                               Condition = parsedCondition,
#                               BioReplicate = myGlobal_condition[,3],
#                               Run = myGlobal_condition[,2],
#                               IsotopeLabelType = rep("L",length(myGlobal_condition[,2])),
#                               stringsAsFactors = F)

# Do edit global annotation to get proper 2 group comparison
# fix(Full_annotation_protein)

# this is only relevant after adapting the annotation file (using fix(annotation))
#write.table(x = Full_annotation_protein, file = "protein_TwoGroup_annotationused.txt")



# same 2grp comparison
(myRef <- as.character(unique(myAnno_protein$Condition)[3]))
quantResTotalProt <- get2GroupAnalysisFromExistingAnnotationFileNoRender(proteinMat = protein,
                                                                                           nrPep = nrPeptides,
                                                                                           mxNAs = maxNAsForComparison,
                                                                                           lg2FCthreshold = absLg2FCThreshold,
                                                                                           qvThreshold = qvalueThreshold,
                                                                                           pathNFileForAnnotation = myAnno_proteinFile,
                                                                                           expName = "totalprot_KOvsCtrl",
                                                                                           resultDirectory = "totalProtein",
                                                                                           DenomReference = myRef,
                                                                                           Pairing = FALSE)






#################################################
#       Merge the GlobalLFQ and Phospho results and Do Candidate calling, NtoCplots, exportTables, MotifPlots
#################################################

# render report first in order to get pseudoData
comparisonName <- "Phospho_twoGroup"
relevantTwoGroupOutput <- quantResPhosphoPep$grp2
pdfName <- paste("pXXXX_", comparisonName, ".pdf", sep="")
rmarkdown::render("Grp2Analysis_NewVersion_phosphopeptide_2019_wMissing.Rmd", params=list(grp=relevantTwoGroupOutput), bookdown::pdf_document2())
file.copy(from = "Grp2Analysis_NewVersion_phosphopeptide_2019_wMissing.pdf", to = pdfName, overwrite=TRUE)
file.remove("Grp2Analysis_NewVersion_phosphopeptide_2019_wMissing.pdf")

#how about there is no pseudo
if (exists("pseudoData")) {
  ActualComplete_phosphoPepQuantResults <- pseudoData
} else {
  ActualComplete_phosphoPepQuantResults <- quantResPhosphoPep$myResultTable
}

# Two Group Comparison
res__merge_allPhosPep <- getDFwithMergedQuantResultsWithLargeMatrix(theGoodMatrix = myMat_good, phosPepsQuantResults = ActualComplete_phosphoPepQuantResults)

# relaxed thresholds here to get some candidates although reduced input data
#my_candidates <- getCandidatesFromMergedMatrix(mergedResultMatrix = res__merge_allPhosPep, my_pModThreshold = pModThreshold, my_absLog2FCThreshold = absLg2FCThreshold)
my_candidates <- getCandidatesFromMergedMatrix(mergedResultMatrix = res__merge_allPhosPep, my_pModThreshold = pModThreshold, my_absLog2FCThreshold = absLg2FCThreshold)
nrow(my_candidates)

# this needs revision for new SRM package till here: 2019-08-20 17:52 -> done! (use _2019 functions)
# relaxed thresholds here to get some candidates although reduced input data
generateTxtFilesForDownstreamAnalysisAndReturnCandidateMatrix_2019(comparisonName = "Comparisond_KOvsCtrl", mergedResultMatrix = res__merge_allPhosPep, my_pModThreshold = pModThreshold, my_absLog2FCThreshold = absLg2FCThreshold, isUniprotDB = FALSE)

my_combo_globalNphosphopep <- combineAndReturnPhosphoMatrix(mergedResultMatrix = res__merge_allPhosPep, globalProteinQuantResults = quantResTotalProt$myResultTable)
generateNtoCProteinPDFsWithPhosphoPeptides_2019(globalNphosphoCombinedNResultMatrix  = my_combo_globalNphosphopep, candidateMatrix = my_candidates, expName = "pXXX_NtoCplots")
write.table(my_combo_globalNphosphopep, file= "pXXX_my_combo_globalNPhospho.txt", row.names = FALSE, quote=FALSE, sep="\t")



# only sequences of the same length go in here
pdf("SequenceLogos_Serine.pdf",8,4)
doSequenceLogoPlotsFromCandidateMatrixForSinglePhosPeptides(merged_phosphoPeptideResultMatrix = my_candidates, phosAcceptorAA = "S", myLogoPlotTitle = "Significant Sequence logos - Serines")
dev.off()


pdf("SequenceLogos_Threonine.pdf",8,4)
doSequenceLogoPlotsFromCandidateMatrixForSinglePhosPeptides(merged_phosphoPeptideResultMatrix = my_candidates, phosAcceptorAA = "T", myLogoPlotTitle = "Significant Sequence logos - Threonine")
dev.off()

