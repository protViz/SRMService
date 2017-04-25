##
##

maxquanttxtdirectory <- ""
evidence <- read.table("evidence.txt")

#
msmsName<- "maxquant/msms.txt"
msms_d <- read.table(msmsName, header=T, sep="\t")

ff <- "maxquant/summary.txt"
summ <- read.table(ff, header=F, sep="\t")
f <- "maxquant/evidence.txt"
evi_d <- read.table(f, header=T, sep="\t")

f <- "maxquant/proteinGroups.txt"
Fulldat <- read.table(f, header=T, sep="\t")
bool_moreThanOnePeptide <- Fulldat$Razor...unique.peptides > 1
f <- "proteinGroups_FGCZ2grp_Intensity.txt"
dat <- read.table(f, header=T, sep="\t",row.names=1)


# will not run..
Stangle('Paired_2grpAnalysis_forProperOrganizedQMatrix_justTheMatrix.Rnw')
Sweave('Paired_2grpAnalysis_forProperOrganizedQMatrix_justTheMatrix.Rnw')
