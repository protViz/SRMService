library(readxl)
library(PECA)
library(plyr)
library(foreach)

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

i <- 1
group1 <- colnames(data)[pairs[1,i]:(pairs[1,i]+2)]
group2 <- colnames(data)[pairs[2,i]:(pairs[2,i]+2)]
group1
group2
head(data)

peca.out.modt <- PECA_df(data, "Protein", group1, group2, test="modt", type="median", normalize=FALSE, paired=FALSE)
peca.out.rots <- PECA_df(data, "Protein", group1, group2, test="rots", type="median", normalize=FALSE, paired=FALSE)


results.ropeca <- foreach(i=1:ncol(pairs), .combine=rbind.data.frame) %do% {
  print(i)
  group1 <- colnames(data)[pairs[1,i]:(pairs[1,i]+2)]
  group2 <- colnames(data)[pairs[2,i]:(pairs[2,i]+2)]

    peca.out <- PECA_df(data, "Protein", group1, group2, test="rots", type="median", normalize=FALSE, paired=FALSE)
  peca.out$Protein <- rownames(peca.out)
  peca.out
}
