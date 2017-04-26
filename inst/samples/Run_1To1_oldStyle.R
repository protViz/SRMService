##
##
library(tools)
maxquanttxtdirectory <- ""
#evidence <- read.table("inst/samples/maxquant_txt/MSQC1/evidence.txt")

#
msmsName<- "inst/samples/maxquant_txt/MSQC1/msms.txt"
msms_d <- read.table(msmsName, header=T, sep="\t")

ff <- "inst/samples/maxquant_txt/MSQC1/summary.txt"
summ <- read.table(ff, header=F, sep="\t")
f <- "inst/samples/maxquant_txt/MSQC1/evidence.txt"
evi_d <- read.table(f, header=T, sep="\t")

f <- "inst/samples/maxquant_txt/MSQC1/proteinGroups.txt"
Fulldat <- read.table(f, header=T, sep="\t")
bool_moreThanOnePeptide <- Fulldat$Razor...unique.peptides > 1

# missing:
# parameters
ff <- "inst/samples/maxquant_txt/MSQC1/parameters.txt"
params <- read.table(ff, header=T, sep="\t")
# peptides
ff <- "inst/samples/maxquant_txt/MSQC1/peptides.txt"
pepts <- read.table(ff, header=T, sep="\t")

# this is for witold to be replaced
f <- "inst/samples/maxquant_txt/MSQC1/proteinGroups_FGCZ2grp_Intensity.txt"
dat <- read.table(f, header=T, sep="\t",row.names=1)

# get it committed

# get Sample QC running
Stangle('inst/samples/MQ_sampleQC_overview.Rnw')
Sweave('inst/samples/MQ_sampleQC_overview.Rnw')
texi2dvi('MQ_sampleQC_overview.tex', pdf = TRUE)

