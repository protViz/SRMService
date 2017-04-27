#R
# January 2014
# some changes along 2015
# around Feb 2015Hubis smart test
# August 2015 .. more options for some functions
# Nov 2015 .. paired t-test
# Dez 2015 .. Multi-group analysis
# Jan 2016 .. Boxplot for 2 groups


# Jonas Grossmann <jg@fgcz.ethz.ch>

###################################################
### Function for ProteinQuantMatrix
###################################################

#library(affy)
library(missForest)
#library(gplots)
#library(limma)
#library(genefilter)
#library(beeswarm)
library(quantable)


#RobustTtest
#CHECKED
RobustTtest.returnpValue <- function(grpOne, grpTwo) {
  obj<-try(t.test(grpOne, grpTwo), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#Add more groups for comparisons
RobustFtest.returnpValue <- function(quantMatrix, groupFactorScenario) {
  obj<-try(rowFtests(quantMatrix, groupFactorScenario), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


# 2grps
Do2grpTtestOnMatrixAndBHcorrReturnAllInternalTrafo = function(ProtQuantMatrix_rn, bool_TrafoHere = TRUE){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]

  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  #ttest
  pValueVector <- vector()
  if (bool_TrafoHere == TRUE) {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustTtest.returnpValue(asinh(grp_1[i,]),asinh(grp_2[i,]))
    }
  } else {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
    }
  }

  #multipleTestingCorrection
  fdrValues <- vector()
  fdrValues <- p.adjust(pValueVector, method="fdr")
  #FC
  log2FCvector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    log2FCvector[i] <- log(grp1_means[i]/grp2_means[i],2)
  }
  proteinNamesPValFCs<- data.frame(matrix(ncol = 4, nrow(ProtQuantMatrix_rn)))
  colnames(proteinNamesPValFCs) <- c("ProtName","pValue","fdr","log2FC")
  proteinNamesPValFCs <- cbind(row.names(ProtQuantMatrix_rn), pValueVector, fdrValues, log2FCvector)
  return(proteinNamesPValFCs)
}





#ProtQuantMatrix_rn <- i_dat
# CHECKED
#ProtQuantMatrix_rn <- mat
#SignificanceThreshold=0.01
#LinFoldChangeThreshold=2
#bool_TrafoHere = TRUE
Do2grpTtestRobustOnMatrixAndBHcorrWithThresholdAndFoldChangeAndReturnOnlySignificantsInternalTrafo = function(ProtQuantMatrix_rn, SignificanceThreshold=0.01, LinFoldChangeThreshold=2, bool_TrafoHere=TRUE){
  FoldChangeThreshold <- abs(log(LinFoldChangeThreshold,2))
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]

  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  #ttest
  pValueVector <- vector()
  if (bool_TrafoHere == TRUE) {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustTtest.returnpValue(asinh(grp_1[i,]),asinh(grp_2[i,]))
    }
  } else {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
    }
  }

  #multipleTestingCorrection
  fdrValues <- vector()
  fdrValues <- p.adjust(pValueVector, method="fdr")
  #FC
  log2FCvector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    log2FCvector[i] <- NA
    #here a try..
    log2FCvector[i] <- log(grp1_means[i]/grp2_means[i],2)
  }
  boolIncludeIt_vec<- vector()
  #df with 3 cols and rownames
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    if (is.finite(log2FCvector[i]) & abs(log2FCvector[i]) > FoldChangeThreshold & is.finite(pValueVector[i]) & pValueVector[i] <  SignificanceThreshold) {
      boolIncludeIt_vec[i] <- TRUE
    } else {
      boolIncludeIt_vec[i] <- FALSE
    }
  }
  proteinNamesPValFCs<- data.frame(matrix(ncol = 4, nrow = sum(boolIncludeIt_vec)))
  proteinNamesPValFCs <- cbind(row.names(ProtQuantMatrix_rn)[boolIncludeIt_vec], pValueVector[boolIncludeIt_vec], fdrValues[boolIncludeIt_vec], log2FCvector[boolIncludeIt_vec])
  return(proteinNamesPValFCs)
}

#ProtQuantMatrix_rn <- mat
# CHECKED
DoCorrelationOnMatrix = function(ProtQuantMatrix_rn, ...){
  PQmat <- na.omit(as.matrix(ProtQuantMatrix_rn))
  ProtcorrMatrix <- cor(PQmat)
  heatmap.2(as.matrix(ProtcorrMatrix),margin=c(10,10),trace="none", ...)
}

#CHECKED
# remove rows with Zeros
RemoveRowsWithZerosFromProtMatrix = function(ProtQuantMatrix_rn){
  #Only Look at Not-Zero-Proteins
  boolZeroProtein <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    if (sum(ProtQuantMatrix_rn[i,] == 0) == 0) {
      boolZeroProtein[i]<-TRUE
    } else {
      boolZeroProtein[i]<-FALSE
    }
  }
  ProtQuantMatrix_rn_noZeros <- ProtQuantMatrix_rn[boolZeroProtein,]
  return(ProtQuantMatrix_rn_noZeros)
}

#CHECKED
# Impute Values for Zeros
ImputeValuesInProtMatrixForRowsWithZeros = function(ProtQuantMatrix_rn){
  NumSamples <- ncol(ProtQuantMatrix_rn)
  #replace zeros with NA
  ProtQuantMatrix_rn[ProtQuantMatrix_rn==0] <- NA
  #  library(missForest)
  dataimpute <- missForest(ProtQuantMatrix_rn)$ximp
  ProtQuantMatrix_rn_imputed<-dataimpute
  colnames(ProtQuantMatrix_rn_imputed) = colnames(ProtQuantMatrix_rn)
  return(ProtQuantMatrix_rn_imputed)
}


# Normalization
# CHECKED
NormalizeWithMedianPQMatrix = function(ProtQuantMatrix_rn){
  maxMedian <- max(apply(na.omit(ProtQuantMatrix_rn[,1:ncol(ProtQuantMatrix_rn)]),2,median))
  ScaleFactorsMed <- 1/(apply(na.omit(ProtQuantMatrix_rn[,1:ncol(ProtQuantMatrix_rn)]),2,median)/maxMedian)
  write.table(ScaleFactorsMed, "appliedScaleFactors_Median.txt",row.names=TRUE,col.names=FALSE)
  nPQmatrix <- ProtQuantMatrix_rn
  for (i in 1:ncol(ProtQuantMatrix_rn)) {
    nPQmatrix[,i]<-ProtQuantMatrix_rn[,i]*ScaleFactorsMed[i]
  }
  return(nPQmatrix)
}






