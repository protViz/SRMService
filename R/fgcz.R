#R
#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following URL
# www.github.com/protViz/SRMService
#

#' Title
#'
#' @param maxquanttxtdirectory
#' @param reportFileBaseName
#' @author W.E. Wolski, J. Grossmann, C. Panse
#' @seealso www.github.com/protViz/SRMService
#' @return nothing
.fgcz_perform_rendering <- function(maxquanttxtdirectory = '.', reportFileBaseName = './MQ_sampleQC_overview'){

   RMD_QC1To1_Old(maxquanttxtdirectory)

   msmsName <- "msms.txt"
   summary <- "summary.txt"
   evidence <- "evidence.txt"
   proteinGroups <- "proteinGroups.txt"
   parameters <- "parameters.txt"
   peptides <- "peptides.txt"

   RnwFile <- paste(reportFileBaseName, "Rnw", sep='.')
   texFile <- paste(reportFileBaseName, "tex", sep='.')

   msms_d <<- read.table(msmsName, header = T, sep="\t")
   summ <<- read.table(summary, header = F, sep="\t")
   evi_d <<- read.table(evidence, header = T, sep="\t")
   Fulldat <<- read.table(proteinGroups, header=T, sep="\t")

   dat <<- Fulldat[,grep("^Intensity\\.", colnames(Fulldat))]
   rownames(dat) <<- Fulldat$Majority.protein.IDs

   bool_moreThanOnePeptide <<- Fulldat$Razor...unique.peptides > 1

   params <<- read.table(parameters, header = TRUE, sep = "\t")

   # peptides
   pepts <<- read.table(peptides, header = TRUE, sep = "\t")

   # this is for witold to be replaced
   #fixedProteingroups <- "proteinGroups_FGCZ2grp_Intensity.txt"
   #dat <- read.table(fixedProteingroups, header=T, sep="\t",row.names=1)

   # get it committed

   # get Sample QC running
   Stangle(RnwFile)
   Sweave(RnwFile)
   tools::texi2dvi(texFIle, pdf = TRUE)
}
