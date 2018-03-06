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

msmsName<- "msms.txt"
summary <- "summary.txt"
evidence <- "evidence.txt"
#proteinGroups <- "proteinGroups.txt"
proteinGroups <- "proteinGroups_withTicksFaked.txt"
parameters <- "parameters.txt"


msms_d <- read.table(msmsName, header=T, sep="\t")
summ <- read.table(summary, header=F, sep="\t")
evi_d <- read.table(evidence, header=T, sep="\t")
Fulldat <- read.table(proteinGroups, header=T, sep="\t")

dat <- Fulldat[,grep("^Intensity\\.", colnames(Fulldat))]
rownames(dat) <- Fulldat$Majority.protein.IDs


bool_moreThanOnePeptide <- Fulldat$Razor...unique.peptides > 1

params <- read.table(parameters, header=T, sep="\t")

# peptides
peptides <- "peptides.txt"
pepts <- read.table(peptides, header=T, sep="\t")

# this is for witold to be replaced
#fixedProteingroups <- "proteinGroups_FGCZ2grp_Intensity.txt"
#dat <- read.table(fixedProteingroups, header=T, sep="\t",row.names=1)



# get Sample QC running
Stangle('MQ_sampleQC_overview.Rnw')
Sweave('MQ_sampleQC_overview.Rnw')
texi2dvi('MQ_sampleQC_overview.tex', pdf = TRUE)

