##
##
library(tools)
maxquanttxtdirectory <- ""
#evidence <- read.table("inst/samples/maxquant_txt/MSQC1/evidence.txt")

#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following URL
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
#

reportFileBaseName <- "fgcz_MQ_QC_report"
# get these from executing shell script
# 2018-08-17 -> convention and cleaning with CP
# e.g. mq_msms_filename
mq_msms_filename <- "msms.txt"
mq_summary_filename <- "summary.txt"
mq_evidence_filename <- "evidence.txt"
mq_proteinGroups_filename <- "proteinGroups.txt"
mq_parameters_filename <- "parameters.txt"
mq_peptides_filename <- "peptides.txt"

texFile <- paste(reportFileBaseName, "tex", sep='.')
RnwFile <- paste(reportFileBaseName, "Rnw", sep=".")

projectID <<- "xxxx"
orderID <<- "xxxx"


# clean names from Intensity and truncate date for better readability
tmp <- tempdir()
rnwFile <- system.file(file.path("OneToOneAnalysis", RnwFile), package="SRMService")
pdfFile <- system.file(file.path("OneToOneAnalysis/images",'LFQ_QC_workflow.pdf'), package="SRMService")

file.copy(rnwFile, tmp, overwrite = TRUE)
file.copy(pdfFile, tmp, overwrite = TRUE)

renderRNW <- file.path(tmp, RnwFile)

# get Sample QC running
wd<-getwd()
setwd(tmp)

Stangle(RnwFile)
Sweave(RnwFile)
texi2dvi(texFile, pdf = TRUE)
setwd(wd)


