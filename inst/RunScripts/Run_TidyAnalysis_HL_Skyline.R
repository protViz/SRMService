rm(list=ls())
library(tidyverse)
library(rlang)
library(conflicted)
library(SRMService)


outdir <- tempdir()

#allDataM <- read.csv(file="d:\\projects\\p2342_JonasZaugg\\order4041\\output/allData.txt", stringsAsFactors = FALSE)
#skylineSRM_HL_data <- allDataM
#usethis::use_data(skylineSRM_HL_data)

allDataM <- SRMService::skylineSRM_HL_data
head(allDataM)

#source("c:/Users/wolski/prog/SRMService/R/tidyMS_R6_AnalysisConfiguration.R")

skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label", qValue="annotation_QValue")
skylineconfig$table$factors[["treatment_c"]] <- "Condition2"
skylineconfig$table$factors[["time_c"]] <- "time"
skylineconfig$parameter$workingIntensityTransform = ""

resData <- setup_analysis(allDataM, skylineconfig)
#sample_analysis_HL <- resData
#usethis::use_data(sample_analysis_HL, overwrite = TRUE)

resData$Area[resData$Area == 0] <- NA

tt <- R6extractValues(skylineconfig)
yaml::write_yaml(tt, file=file.path(outdir,"testSkyline.yml"))

proteinIDsymbol <- sym(skylineconfig$table$hierarchyKeys()[1])
xnested <- resData %>% group_by(UQ(proteinIDsymbol)) %>% nest()

linePlotHierarchy_configuration(xnested$data[[2]], xnested$protein_Id[[2]], skylineconfig)
linePlotHierarchy_configuration(xnested$data[[3]], xnested$protein_Id[[3]], skylineconfig, separate = TRUE)

figs <- xnested %>%
  mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig, separate = TRUE))

if(1){
  pdf(file.path(outdir,"allProteins.pdf"), width = 10, height = 10)
  invisible(lapply(figs$plot[1:3], print))
  dev.off()
}

#rmarkdown::render("SRMReport.Rmd", params = list(data=resData, configuration = skylineconfig), envir = new.env(), output_file = "SRMReport.pdf")
HLData <- spreadValueVarsIsotopeLabel(resData, skylineconfig)
HLData <- HLData %>% mutate(log2L_log2H = log2(light_Area)- log2(heavy_Area))
HLData$Isotope.Label <- "L/H"

skylineconfigHL <- skylineconfig
skylineconfigHL$table$workIntensity = "log2L_log2H"
skylineconfigHL$parameter$workingIntensityTransform = "log"
skylineconfigHL$table$isotopeLabel

xnested <- HLData %>% group_by(UQ(proteinIDsymbol)) %>% nest()
linePlotHierarchy_configuration(xnested$data[[2]], xnested$protein_Id[[2]], skylineconfigHL)

HLfigs <- xnested %>% mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig))

if(1){
  pdf(file.path(outdir,"allHLProteins.pdf"), width = 10, height = 10)
  invisible(lapply(HLfigs$plot[1:3], print))
  dev.off()
}


options(warn=0)
HLfigs2 <- HLfigs %>% mutate(spreadMatrix = map(data, extractIntensities, skylineconfig))
HLfigs2 <- HLfigs2 %>% mutate(medpolishRes = map(spreadMatrix, medpolishPly))

#reestablishCondition(HLfigs2$data[[2]], HLfigs2$medpolishRes[[2]], skylineconfig )

HLfigs3 <- HLfigs2 %>% mutate(medpolishRes = map2(data,medpolishRes,reestablishCondition , skylineconfig))

HLfigs3 <- HLfigs3 %>% mutate(figsMed = map2(plot, medpolishRes, linePlotHierarchy_QuantLine, "medpolish", skylineconfig))


pdf(file.path(outdir,"allProteinsWithMed.pdf"), width = 10, height = 10)
invisible(lapply(HLfigs3$figsMed[2:5], print))
dev.off()

prots <- HLfigs3 %>% dplyr::select(as.character(proteinIDsymbol),  medpolishRes) %>% unnest()

write_tsv(prots, path = file.path(outdir,"proteins.tsv"))

pdf(file.path(outdir,"vsTimeYfixed.pdf"), width = 10, height = 10)
ggplot(prots, aes(x=time_c, y=medpolish, group=treatment_c, color=treatment_c )) +
  geom_point() +
  geom_line() +
  facet_wrap(~protein_Id) +
  theme_classic()
dev.off()


pdf(file.path(outdir,"vsTimeYfree.pdf"), width = 10, height = 10)
ggplot(prots, aes(x=time_c, y=medpolish, group=treatment_c, color=treatment_c )) +
  geom_point() +
  geom_line() +
  facet_wrap(~protein_Id, scales="free_y" ) +
  theme_classic()
dev.off()





