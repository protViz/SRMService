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

# get these from executing shell script
projectNumber <- "p1000"
OID <- "OID1234"


msmsName<- "msms.txt"
summary <- "summary.txt"
evidence <- "evidence.txt"
proteinGroups <- "proteinGroups.txt"
parameters <- "parameters.txt"
peptides <- "peptides.txt"

# 
msmsName <- system.file(file.path("samples/maxquant_txt/MSQC1",msmsName),package = "SRMService")
summary <- system.file(file.path("samples/maxquant_txt/MSQC1",summary),package = "SRMService")
evidence <- system.file(file.path("samples/maxquant_txt/MSQC1",evidence),package = "SRMService")
proteinGroups <- system.file(file.path("samples/maxquant_txt/MSQC1",proteinGroups),package = "SRMService")
parameters <- system.file(file.path("samples/maxquant_txt/MSQC1",parameters),package = "SRMService")
peptides <- system.file(file.path("samples/maxquant_txt/MSQC1",peptides),package = "SRMService")

# msmsName <- system.file(file.path("samples/maxquant_txt/maxquant",msmsName),package = "SRMService")
# summary <- system.file(file.path("samples/maxquant_txt/maxquant",summary),package = "SRMService")
# evidence <- system.file(file.path("samples/maxquant_txt/maxquant",evidence),package = "SRMService")
# proteinGroups <- system.file(file.path("samples/maxquant_txt/maxquant",proteinGroups),package = "SRMService")
# parameters <- system.file(file.path("samples/maxquant_txt/maxquant",parameters),package = "SRMService")
# peptides <- system.file(file.path("samples/maxquant_txt/maxquant",peptides),package = "SRMService")

msms_d <- read.table(msmsName, header=T, sep="\t")
summ <- read.table(summary, header=F, sep="\t")
evi_d <- read.table(evidence, header=T, sep="\t")
params <- read.table(parameters, header=T, sep="\t")
pepts <- read.table(peptides, header=T, sep="\t")

Fulldat <- read.table(proteinGroups, header=T, sep="\t")
dat <- Fulldat[,grep("^Intensity\\.", colnames(Fulldat))]
rownames(dat) <- Fulldat$Majority.protein.IDs

# clean names from Intensity and truncate date for better readability
noIntensityNames <- gsub(pattern = "Intensity.", replacement = "", x = colnames(dat))
stoopidNames <- gsub(pattern = "^m", replacement = "", x = noIntensityNames)
noDateNames <- gsub(pattern = "^[[:digit:]]+_", replacement = "", x = stoopidNames)
colnames(dat) <- noDateNames

# this is for witold to be replaced
#fixedProteingroups <- "proteinGroups_FGCZ2grp_Intensity.txt"
#dat <- read.table(fixedProteingroups, header=T, sep="\t",row.names=1)

tmp <- tempdir()
rnwFile <- system.file(file.path("OneToOneAnalysis",'MQ_sampleQC_overview.Rnw'), package="SRMService")
pdfFile <- system.file(file.path("OneToOneAnalysis/images",'LFQ_QC_workflow.pdf'), package="SRMService")

file.copy(rnwFile, tmp, overwrite = TRUE)
file.copy(pdfFile, tmp, overwrite = TRUE)

renderRNW <- file.path(tmp, 'MQ_sampleQC_overview.Rnw')

# get Sample QC running
wd<-getwd()
setwd(tmp)

Stangle('MQ_sampleQC_overview.Rnw')
Sweave('MQ_sampleQC_overview.Rnw')
texi2dvi('MQ_sampleQC_overview.tex', pdf = TRUE)
setwd(wd)


