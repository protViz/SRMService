rm(list=ls())
library(tidyverse)
library(readxl)
library(rlang)
library(yaml)

data = 2
if(data == 0){
  annotation <- read_excel(path="3DLiver_Isoniazid_Annotation.xlsx")
  colnames(annotation) <- make.names(colnames(annotation))
  skylineDataX <- read.csv("QuantificationPRM_3DLiver_Isoniazid_20180320_NSK.csv", stringsAsFactors = FALSE)
  outdir = "outdir_3DLiver_Isoniazid"
}else if(data == 1){
  annotation <- read_excel(path="3DLiver_Cyclosporin_Annotation.xlsx")
  colnames(annotation) <- make.names(colnames(annotation))
  skylineDataX <- read.csv("QuantificationPRM_3DLiver_Cyclosporin_20180417_NSK.csv", stringsAsFactors = FALSE)
  outdir = "outdir_3DLiver_Cyclosporin"
}else if(data == 2){
  annotation <- read_excel(path="3DLiver_ValproicAcid_Annotation.xlsx")
  colnames(annotation) <- make.names(colnames(annotation))
  skylineDataX <- read.csv("QuantificationPRM_3DLiver_ValproicAcid_20180618_NSK.csv", stringsAsFactors = FALSE)
  outdir = "outdir_3DLiver_ValproicAcid"
}

if(!dir.exists(outdir)){
  dir.create(outdir)
}

skylineData <- inner_join(skylineDataX, annotation, by = c("Replicate.Name"= "Raw.file_PRM"))

length(unique(annotation$Raw.file_PRM))
length(unique(skylineDataX$Replicate.Name))
length(unique(skylineData$Replicate.Name))
skylineData %>% select(Sampling.Time.Point, Replicate.Name) %>% distinct() %>% arrange(Sampling.Time.Point)

source("c:/Users/wolski/prog/SRMService/R/tidyMS_R6_AnalysisConfiguration.R")

skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label.Type", qValue="Detection.Q.Value")
skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"

table <- skylineconfig$table

tt <- R6extractValues(skylineconfig)
write_yaml(tt, file="testSkyline.yml")

library(usethis)

saveRDS(skylineData,file = "c:/Users/wolski/prog/SRMService/inst/skylineData.rda")
head(skylineData)
resData <- setupDataFrame(skylineData, skylineconfig)

xx <- unique(resData$Time)
xxord <- order(as.numeric(gsub("T","",xx)))

resData$Time <- parse_factor(resData$Time, unique(resData$Time)[xxord])

resData$Area[resData$Area == 0] <- NA



### Generate overview plots
proteinIDsymbol <- sym(names(skylineconfig$table$hierarchy[1]))
xnested <- resData %>% group_by(UQ(proteinIDsymbol)) %>% nest()
figs <- xnested %>% mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig))

if(1){
  pdf(file.path(outdir,"allProteinsTransitions.pdf"), width = 10, height = 10)
  invisible(lapply(figs$plot, print))
  dev.off()
}

rmarkdown::render("SRMReport.Rmd", output_format = "html_document", params = list(data = resData, configuration = skylineconfig),envir = new.env())
file.copy("SRMReport.html",outdir )
file.remove("SRMReport.html")
rmarkdown::render("SRMReport.Rmd", output_format = "pdf_document", params = list(data = resData, configuration = skylineconfig),envir = new.env())
file.copy("SRMReport.pdf",outdir )
file.remove("SRMReport.pdf")

resDataLog <- resData %>% mutate(log2Area = log2(Area))
skylineconfig$table$workIntensity = "log2Area"
skylineconfig$parameter$workingIntensityTransform = "log"

xnested <- resDataLog %>% group_by_at(names(table$hierarchy)[1]) %>% nest()
figs <- xnested %>% mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig))
figs2 <- figs %>% mutate(spreadMatrix = map(data, extractIntensities, skylineconfig))
figs2 <- figs2 %>% mutate(medpolishRes = map(spreadMatrix, medpolishPly))
figs3 <- figs2 %>% mutate(medpolishRes = map2(data,medpolishRes,reestablishCondition, skylineconfig ))


linePlotHierarchy_QuantLine(figs3$plot[[1]], figs3$medpolishRes[[1]], "medpolish", skylineconfig)
figs3 <- figs3 %>% mutate(figsMed = map2(plot, medpolishRes, linePlotHierarchy_QuantLine, "medpolish" , skylineconfig ))
head(figs3)

pdf(file.path( outdir , "allProteinsWithMed.pdf"), width = 10, height = 10)
invisible(lapply(figs3$figsMed, print))
dev.off()

protIntensity <- figs3 %>% select(names(table$hierarchy)[1], medpolishRes) %>% unnest()
CiRT <- protIntensity %>% filter(protein_Id == "CiRT standards")
dim(protIntensity)

proteinIntensity <- protIntensity %>%
  inner_join(CiRT, by= setdiff(names(protIntensity), c("protein_Id","medpolish")), suffix = c("",".CiRT")) %>%
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
#protRez3 %>% select( Protein.Name, ba) %>% unnest()
#broom::tidy(protRez3$anova)

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
invisible(lapply(protRez3$plot, print))
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
invisible(lapply(protRez3$plotMedian, print))
dev.off()

