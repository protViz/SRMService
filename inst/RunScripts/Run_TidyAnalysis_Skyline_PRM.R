rm(list=ls())
library(tidyverse)
library(readxl)
library(rlang)
library(yaml)
library(conflicted)

outdir <- tempdir()
data(skylinePRMSampleData)

skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label.Type", qValue="Detection.Q.Value")
skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#usethis::use_data(skylineconfig, overwrite = TRUE)
resData <- setup_analysis(skylinePRMSampleData, skylineconfig)

xx <- unique(resData$Time)
xxord <- order(as.numeric(gsub("T","",xx)))

resData$Time <- parse_factor(resData$Time, unique(resData$Time)[xxord])
resData$Area[resData$Area == 0] <- NA



### Generate overview plots
proteinIDsymbol <- sym(names(skylineconfig$table$hierarchy)[1])
xnested <- resData %>% group_by(UQ(proteinIDsymbol)) %>% nest()
figs <- xnested %>% mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig))

print(figs$plot[[1]])

if(1){
  pdf(file.path(outdir,"allProteinsTransitions.pdf"), width = 10, height = 10)
  invisible(lapply(figs$plot[1:3], print))
  dev.off()
}

resDataLog <- resData %>% mutate(log2Area = log2(Area))
skylineconfig$table$workIntensity = "log2Area"
skylineconfig$parameter$workingIntensityTransform = "log"

#source("c:/Users/wolski/prog/SRMService/R/tidyMS_R6_AnalysisConfiguration.R")

xnested <- resDataLog %>% group_by(UQ(proteinIDsymbol)) %>% nest()
figs <- xnested %>% mutate(plotlog = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig))
print(figs$plotlog[[1]])

figs2 <- figs %>% mutate(spreadMatrix = map(data, extractIntensities, skylineconfig))
figs2 <- figs2 %>% mutate(medpolishRes = map(spreadMatrix, medpolishPly))
figs3 <- figs2 %>% mutate(medpolishRes = map2(data,medpolishRes,reestablishCondition, skylineconfig ))


linePlotHierarchy_QuantLine(figs3$plotlog[[1]], figs3$medpolishRes[[1]], "medpolish", skylineconfig)

figs3 <- figs3 %>% mutate(figsMed = map2(plotlog, medpolishRes, linePlotHierarchy_QuantLine, "medpolish" , skylineconfig ))
head(figs3)
print(figs3$figsMed[[1]])

pdf(file.path( outdir , "allProteinsWithMed.pdf"), width = 10, height = 10)
invisible(lapply(figs3$figsMed[1:3], print))
dev.off()

table <- skylineconfig$table
protIntensity <- figs3 %>% select(names(table$hierarchy)[1], medpolishRes) %>% unnest()
CiRT <- protIntensity %>% dplyr::filter(protein_Id == "CiRT standards")
dim(protIntensity)

proteinIntensity <- protIntensity %>%
  inner_join(CiRT, by= dplyr::setdiff(names(protIntensity), c("protein_Id","medpolish")), suffix = c("",".CiRT")) %>%
  mutate(log2Med_log2MedCiRT = medpolish - medpolish.CiRT)


p <- ggplot(proteinIntensity, aes(x =sampleName , y=log2Med_log2MedCiRT, group=protein_Id, color=protein_Id)) +
  geom_line() +
  facet_grid(~Time, scales = "free_x") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  ylab(expression(log[2](frac(P,CiRT))))
p

pdf(file.path(outdir,"allProteinsOnePlot.pdf"))
print(p)
dev.off()

write_tsv(proteinIntensity,path=file.path(outdir, "ProteinQuants.tsv"))

protRez3 <- proteinIntensity %>% group_by(protein_Id) %>% nest()
tmp <- function(x){ anova(lm(log2Med_log2MedCiRT ~ Time, data=x)) }

protRez3 <- protRez3 %>% mutate(anova = map( data, tmp))
protRez3 <- protRez3 %>% mutate(ba = map(anova, broom::tidy))

plotProt <- function(data, title){
  ggplot(data, aes(x = Time, y=log2Med_log2MedCiRT )) +
    geom_point(fill="red", color=2, size=2) +
    stat_summary(fun.y=mean, geom="line", aes(group=1), col="blue") +
    stat_summary(fun.y=mean, geom="point", col="blue") + ggtitle(title) +
    theme_classic() + ylab(expression(log[2](frac(P,CiRT))))
}

plotProt(protRez3$data[[3]], protRez3$protein_Id[[3]])
protRez3 <- protRez3 %>% mutate(plot = map2(data,protein_Id, plotProt))

pdf(file.path(outdir,"allProtFigsNorm.pdf"))
invisible(lapply(protRez3$plot[1:3], print))
dev.off()


plotProtMedian <- function(data, title){
  ggplot(data, aes(x = Time, y=log2Med_log2MedCiRT )) +
    geom_point(fill="red", color=2, size=2) +
    stat_summary(fun.y=median, geom="line", aes(group=1), col="blue") +
    stat_summary(fun.y=median, geom="point", col="blue") + ggtitle(title) +
    theme_classic() + ylab(expression(log[2](frac(P,CiRT))))
}


protRez3 <- protRez3 %>% mutate(plotMedian = map2(data,protein_Id, plotProtMedian))
pdf(file.path(outdir,"allProtFigsNormMedian.pdf"))
invisible(lapply(protRez3$plotMedian[1:3], print))
dev.off()

