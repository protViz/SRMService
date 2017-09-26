library(readxl)
library(genefilter)
library(plyr)
library(foreach)
library(MSstats)
library(PECA)

#--------------------------------------------------
# Read and prepare data
#--------------------------------------------------

# Read data
data <- read_excel("Data/mcp.M114.044305-3.xlsx", sheet=2)
colnames(data) <- c("Protein","Peptide","Charge","Sample1_R1","Sample1_R2","Sample1_R3","Sample2_R1","Sample2_R2","Sample2_R3","Sample3_R1","Sample3_R2","Sample3_R3","Sample4_R1","Sample4_R2","Sample4_R3","Sample5_R1","Sample5_R2","Sample5_R3","Sample6_R1","Sample6_R2","Sample6_R3","Sample7_R1","Sample7_R2","Sample7_R3","Sample8_R1","Sample8_R2","Sample8_R3")
data <- data[-1,]
data <- data.frame(data)
head(data)

# Fill in protein identifier and remove empty rows
protein <- "X"
pb <- txtProgressBar(min=1, max=nrow(data), style=3)
for (i in 1:nrow(data)) {
  if (is.na(data[i,"Protein"])) {
    data[i,"Protein"] <- protein
  } else {
    protein <- data[i,"Protein"]
  }
  setTxtProgressBar(pb, i)
}
close(pb)
data <- data[!is.na(data[,"Peptide"]),]

# Convert values to numeric
for (i in 3:ncol(data)) {
  data[,i] <- as.numeric(data[,i])
}



#--------------------------------------------------
# ROPECA
#--------------------------------------------------
pairs <- combn(seq(from=4,to=22,by=3), 2)
results.ropeca <- foreach(i=1:ncol(pairs), .combine=rbind.data.frame) %do% {
  group1 <- colnames(data)[pairs[1,i]:(pairs[1,i]+2)]
  group2 <- colnames(data)[pairs[2,i]:(pairs[2,i]+2)]
  peca.out <- PECA_df(data, "Protein", group1, group2, test="rots", type="median", normalize=FALSE, paired=FALSE)
  peca.out$Protein <- rownames(peca.out)
  peca.out
}

#--------------------------------------------------
# MSstats
#--------------------------------------------------

# Function to get intensitites in MSstats format
getIntensities <- function(condition, bioc, run, sample) {
  RawData <- data.frame(
    ProteinName=data[,"Protein"],
    PeptideSequence=data[,"Peptide"],
    PrecursorCharge=data[,"Charge"],
    FragmentIon=NA,
    ProductCharge=NA,
    IsotopeLabelType="L",
    Condition=condition,
    BioReplicate=bioc,
    Run=run,
    Intensity=data[,sample]
  )
  colnames(RawData) <- c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
  return(RawData)
}

# MSstats for the pairs
pairs <- combn(seq(from=4,to=22,by=3), 2)
results.msstats <- foreach(i=1:ncol(pairs), .combine=rbind.data.frame) %do% {
  data.msstats <- rbind(
    getIntensities(1,"ReplA",1,pairs[1,i]),
    getIntensities(1,"ReplA",2,pairs[1,i]+1),
    getIntensities(1,"ReplA",3,pairs[1,i]+2),
    getIntensities(2,"ReplB",4,pairs[2,i]),
    getIntensities(2,"ReplB",5,pairs[2,i]+1),
    getIntensities(2,"ReplB",6,pairs[2,i]+2)
  )
  QuantData <- dataProcess(data.msstats, normalization=FALSE)
  comparison <- matrix(c(-1,1),nrow=1)
  row.names(comparison) <- "1-2"
  msstats.out <- groupComparison(contrast.matrix=comparison, data=QuantData)
  msstats.out <- data.frame(p=msstats.out$ComparisonResult$pvalue, p.fdr=msstats.out$ComparisonResult$adj.pvalue, protein=msstats.out$ComparisonResult$Protein, stringsAsFactors=FALSE)
  msstats.out
}

#--------------------------------------------------
# t-test
#--------------------------------------------------

# Summarize to proteins
data <- ddply(data, .(Protein), numcolwise(mean,na.rm=TRUE))

# t-test for the pairs
pairs <- combn(seq(from=3,to=21,by=3), 2)
results.ttest <- foreach(i=1:ncol(pairs), .combine=rbind.data.frame) %do% {
  group1 <- pairs[1,i]:(pairs[1,i]+2)
  group2 <- pairs[2,i]:(pairs[2,i]+2)
  labels <- factor(c(1,1,1,2,2,2))
  ttest.out <- rowttests(as.matrix(log2(data[,c(group1,group2)]+1)), fac=labels)
  ttest.out$p.fdr <- p.adjust(ttest.out$p.value, method="BH")
  ttest.out$protein <- data$Protein
  ttest.out
}
