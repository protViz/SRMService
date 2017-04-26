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

library(affy)
library(missForest)
library(gplots)
library(limma)
library(genefilter)
library(beeswarm)
library(quantable)


# 2grp
# CHECKED
Do2grpTtestOnMatrixWithThresholdAndReturnOnlySignificants = function(ProtQuantMatrix_rn, SignificanceThreshold){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
  }
  
  sigThreshold <- SignificanceThreshold
  bool_pVector<- vector()
  outVector<- data.frame(matrix(ncol = 1, nrow = sum(bool_pVector)))
  bool_pVector <- pValueVector<sigThreshold
  row.names(outVector) <- row.names(ProtQuantMatrix_rn)[bool_pVector]
  outVector <- pValueVector[bool_pVector]
  return(outVector)
}


Do2grpTtestOptionallyPairedOnMatrixWithThresholdAndReturnOnlySignificants = function(ProtQuantMatrix_rn, SignificanceThreshold, Bool_pairedOrNot=FALSE){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(grp_1[i,],grp_2[i,], Bool_pairedOrNot)
  }
  
  sigThreshold <- SignificanceThreshold
  bool_pVector<- vector()
  outVector<- data.frame(matrix(ncol = 1, nrow = sum(bool_pVector)))
  bool_pVector <- pValueVector<sigThreshold
  row.names(outVector) <- row.names(ProtQuantMatrix_rn)[bool_pVector]
  outVector <- pValueVector[bool_pVector]
  return(outVector)
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




# 2grps
Do2grpTtestOptionallyPairedOnMatrixAndBHcorrReturnAllInternalTrafo = function(ProtQuantMatrix_rn, bool_TrafoHere = TRUE, Bool_pairedOrNot=FALSE){
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
      pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(asinh(grp_1[i,]),asinh(grp_2[i,]), pairedTRUEorFalse = Bool_pairedOrNot)
    }
  } else {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(grp_1[i,],grp_2[i,], pairedTRUEorFalse = Bool_pairedOrNot)
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



#MultiGroupRobustFtestOnMatrixWBHCorrWithpVThresholdReturnOnlySignificantsAndInternalTrafo(ProtQuantMatrix_rn = n_i_dat, GroupingAsFactors = fGroups,SignificanceThreshold = 0.05,bool_TrafoHere = TRUE)
#ProtQuantMatrix_rn = n_i_dat
#GroupingAsFactors = fGroups
#SignificanceThreshold = 0.05
#bool_TrafoHere = TRUE

DoMultiGroupRobustFtestOnMatrixWBHCorrWithpVThresholdReturnOnlySignificantsAndInternalTrafo = function(ProtQuantMatrix_rn, GroupingAsFactors,  SignificanceThreshold=0.01, bool_TrafoHere = TRUE){
  #F test
  pValueVector <- vector()
  matForTest <- as.matrix(ProtQuantMatrix_rn)
  if (bool_TrafoHere == TRUE) {
    matForTest <- asinh(matForTest)
  }
  pValueVector = RobustFtest.returnpValue(quantMatrix = matForTest, groupFactorScenario = GroupingAsFactors)
  fdrValues <- p.adjust(pValueVector, method="fdr")
  boolIncludeIt_vec <- is.finite(pValueVector) & pValueVector <  SignificanceThreshold
  proteinNamesPValFCs<- data.frame("ProteinGroup"=rownames(ProtQuantMatrix_rn),
                                   "pValue"=pValueVector,
                                   "FDR(adjusted pValue)"=fdrValues,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE)[boolIncludeIt_vec, ,drop=FALSE]
  return(proteinNamesPValFCs)
}

#l <- MultiGroupRobustFtestOnMatrixWBHCorrWithpVThresholdReturnOnlySignificantsAndInternalTrafo(ProtQuantMatrix_rn = n_i_dat, GroupingAsFactors = fGroups,SignificanceThreshold = 0.05,bool_TrafoHere = TRUE)


DoMultiGroupRobustFtestOnMatrixWBHCorrWithpVThresholdReturnFullListInternalTrafo = function(ProtQuantMatrix_rn, GroupingAsFactors, bool_TrafoHere = TRUE){
  #F test
  pValueVector <- vector()
  matForTest <- as.matrix(ProtQuantMatrix_rn)
  if (bool_TrafoHere == TRUE) {
    matForTest <- asinh(matForTest)
  }
  pValueVector = RobustFtest.returnpValue(quantMatrix = matForTest, groupFactorScenario = GroupingAsFactors)
  fdrValues <- p.adjust(pValueVector, method="fdr")
  proteinNamesPValFCs<- data.frame("ProteinGroup"=rownames(ProtQuantMatrix_rn),
                                   "pValue"=pValueVector,
                                   "FDR(adjusted pValue)"=fdrValues,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE)
  return(proteinNamesPValFCs)
}





#i <- 3
#asinh(grp_1[i,]),asinh(grp_2[i,])
#RobustPairedTtest
#pairedTRUEorFalse <- TRUE
#CHECKED
#t.test(as.numeric(asinh(grp_1[i,])), as.numeric(asinh(grp_2[i,])), paired = TRUE)
#t.test(c(0,1,2,3,4), c(1,2,3,4,5), paired = TRUE)
RobustPairedOrUnpairedTtest.returnpValue <- function(grpOne, grpTwo, pairedTRUEorFalse=FALSE) {
  obj<-try(t.test(as.numeric(grpOne), as.numeric(grpTwo), paired = pairedTRUEorFalse), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
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


#LinFoldChangeThreshold=2
#bool_TrafoHere = TRUE
#ProtQuantMatrix_rn <- mat
Bool_pairedOrNot=TRUE
DoOptionallyPairedOrUnpaired2grpTtestRobustOnMatrixAndBHcorrWithThresholdAndFoldChangeAndReturnOnlySignificantsInternalTrafo = function(ProtQuantMatrix_rn, SignificanceThreshold=0.01, LinFoldChangeThreshold=2, bool_TrafoHere = TRUE, Bool_pairedOrNot=FALSE){
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
      pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(asinh(grp_1[i,]),asinh(grp_2[i,]), pairedTRUEorFalse = Bool_pairedOrNot)
    }
  } else {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(grp_1[i,],grp_2[i,], pairedTRUEorFalse = Bool_pairedOrNot)
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



### this is a rewrite of Do2grpTtestRobustOnMatrixAndBHcorrWithThresholdAndFoldChangeAndReturnOnlySignificantsInternalTrafo
## it gets the same results; but does not cast the output into a character matrix with the cbind
HubisSmartTestFor2GrpTtestRobstOnMatrixWBHCorrWithThresholdAndFCAndReturnOnlySignificantsAndInternalTrafo = function(ProtQuantMatrix_rn, SignificanceThreshold=0.01, LinFoldChangeThreshold=2, bool_TrafoHere = TRUE){
  FoldChangeThreshold <- abs(log2(LinFoldChangeThreshold))
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #mean vectors
  grp1_means <- rowMeans(grp_1)
  grp2_means <- rowMeans(grp_2)
  #ttest
  pValueVector <- vector()
  if (bool_TrafoHere == TRUE) {
    grp_1 = asinh(grp_1)
    grp_2 = asinh(grp_2)
  }
  pValueVector = sapply(1:nrow(ProtQuantMatrix_rn), function(i){RobustTtest.returnpValue(grp_1[i,], grp_2[i,])})
  fdrValues <- p.adjust(pValueVector, method="fdr")
  log2FCvector <- log2(grp1_means/grp2_means)
  boolIncludeIt_vec <- is.finite(log2FCvector) & abs(log2FCvector) > FoldChangeThreshold & is.finite(pValueVector) & pValueVector <  SignificanceThreshold
  proteinNamesPValFCs<- data.frame("ProteinGroup"=rownames(ProtQuantMatrix_rn),
                                   "pValue"=pValueVector,
                                   "FDR(adjusted pValue)"=fdrValues,
                                   "log2(mean(grp1)/mean(grp2))"=log2FCvector,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE)[boolIncludeIt_vec, ,drop=FALSE]
  return(proteinNamesPValFCs)
}


### this is a rewrite of Do2grpTtestRobustOnMatrixAndBHcorrWithThresholdAndFoldChangeAndReturnOnlySignificantsInternalTrafo
## it gets the same results; but does not cast the output into a character matrix with the cbind
HubisSmartTestFor2GrpTtestOptionallyPairedRobstOnMatrixWBHCorrWithThresholdAndFCAndReturnOnlySignificantsAndInternalTrafo = function(ProtQuantMatrix_rn, SignificanceThreshold=0.01, LinFoldChangeThreshold=2, bool_TrafoHere = TRUE, Bool_pairedOrNot=FALSE){
  FoldChangeThreshold <- abs(log2(LinFoldChangeThreshold))
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #mean vectors
  grp1_means <- rowMeans(grp_1)
  grp2_means <- rowMeans(grp_2)
  #ttest
  pValueVector <- vector()
  if (bool_TrafoHere == TRUE) {
    grp_1 = asinh(grp_1)
    grp_2 = asinh(grp_2)
  }
  pValueVector = sapply(1:nrow(ProtQuantMatrix_rn), function(i){RobustPairedOrUnpairedTtest.returnpValue(grp_1[i,], grp_2[i,], pairedTRUEorFalse = Bool_pairedOrNot)})
  fdrValues <- p.adjust(pValueVector, method="fdr")
  log2FCvector <- log2(grp1_means/grp2_means)
  boolIncludeIt_vec <- is.finite(log2FCvector) & abs(log2FCvector) > FoldChangeThreshold & is.finite(pValueVector) & pValueVector <  SignificanceThreshold
  proteinNamesPValFCs<- data.frame("ProteinGroup"=rownames(ProtQuantMatrix_rn),
                                   "pValue"=pValueVector,
                                   "FDR(adjusted pValue)"=fdrValues,
                                   "log2(mean(grp1)/mean(grp2))"=log2FCvector,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE)[boolIncludeIt_vec, ,drop=FALSE]
  return(proteinNamesPValFCs)
}




#ProtQuantMatrix_rn <- mat
# CHECKED
DoCorrelationOnMatrix = function(ProtQuantMatrix_rn, ...){
  PQmat <- na.omit(as.matrix(ProtQuantMatrix_rn))
  ProtcorrMatrix <- cor(PQmat)
  heatmap.2(as.matrix(ProtcorrMatrix),margin=c(10,10),trace="none", ...)
}

#
#ProtQuantMatrix_rn <- asinh(n_i_dat)
#groupingFactorsForColoring <- fGroups
#rm(ProtQuantMatrix_rn, groupingFactorsForColoring)

DoCorrelationOnMatrixWsimpleHeatmap = function(ProtQuantMatrix_rn, groupingFactorsForColoring, ...){
  PQmat <- as.matrix(ProtQuantMatrix_rn)
  ProtcorrMatrix <- cor(PQmat)
  simpleheatmap(as.matrix(ProtcorrMatrix), ColSideColors=as.character(as.integer(groupingFactorsForColoring)))
}


# FC calc
DoLog2FCCalc = function(ProtQuantMatrix_rn){
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
  
  log2FCvector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    log2FCvector[i] <- log(grp1_means[i]/grp2_means[i],2) 
  }
  return(log2FCvector)
}


#CHECKED
DoVolcanoPlotWithFCnSigThreshold = function(ProtQuantMatrix_rn, SignificanceThreshold,FoldChangeThreshold, boolTrafo=TRUE){
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
  #check for Trafo
  if (boolTrafo == TRUE) {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustTtest.returnpValue(asinh(grp_1[i,]),asinh(grp_2[i,]))
    }
  } else {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
    }
  }
  log2FCvector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    log2FCvector[i] <- log(grp1_means[i]/grp2_means[i],2) 
  }
  
  SigBoolVec <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    if (is.finite(log2FCvector[i]) & abs(log2FCvector[i]) > log(FoldChangeThreshold,2) & pValueVector[i] <  SignificanceThreshold) {
      SigBoolVec[i]<-TRUE
    } else {
      SigBoolVec[i]<-FALSE
    }
  }
  plot(-1*log(pValueVector,10)~log2FCvector, xlab="log2(grp1/grp2)", ylab="-log10(pValue)",cex=0.5, xlim=c(-6,6), main="Volcano Plot", pch=".")
  points(log2FCvector[SigBoolVec],-1*log(pValueVector[SigBoolVec],10), type="p", col="green")
  abline(v=median(na.omit(log2FCvector)), col="lightblue", lty=2, lwd=2)
  abline(v=0, col="blue")
  abline(v=log(FoldChangeThreshold,2), col="red")
  abline(v=-log(FoldChangeThreshold,2), col="red")
  abline(h=log(SignificanceThreshold,10)*-1, col="green")
  p_names <- row.names(ProtQuantMatrix_rn)
  proteinNamesPValFCs <- cbind(as.character(p_names), pValueVector, log2FCvector)
  return(proteinNamesPValFCs)
}



#CHECKED
DoVolcanoOptionallyPairedPlotWithFCnSigThreshold = function(ProtQuantMatrix_rn, SignificanceThreshold,FoldChangeThreshold, boolTrafo=TRUE, Bool_pairedOrNot=FALSE){
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
  #check for Trafo
  if (boolTrafo == TRUE) {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(asinh(grp_1[i,]),asinh(grp_2[i,]),pairedTRUEorFalse = Bool_pairedOrNot)
    }
  } else {
    for (i in 1:nrow(ProtQuantMatrix_rn)) {
      pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(grp_1[i,],grp_2[i,],pairedTRUEorFalse = Bool_pairedOrNot)
    }
  }
  log2FCvector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    log2FCvector[i] <- log(grp1_means[i]/grp2_means[i],2) 
  }
  
  SigBoolVec <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    if (is.finite(log2FCvector[i]) & abs(log2FCvector[i]) > log(FoldChangeThreshold,2) & pValueVector[i] <  SignificanceThreshold) {
      SigBoolVec[i]<-TRUE
    } else {
      SigBoolVec[i]<-FALSE
    }
  }
  plot(-1*log(pValueVector,10)~log2FCvector, xlab="log2(grp1/grp2)", ylab="-log10(pValue)",cex=0.5, xlim=c(-6,6), main="Volcano Plot", pch=".")
  points(log2FCvector[SigBoolVec],-1*log(pValueVector[SigBoolVec],10), type="p", col="green")
  abline(v=median(na.omit(log2FCvector)), col="lightblue", lty=2, lwd=2)
  abline(v=0, col="blue")
  abline(v=log(FoldChangeThreshold,2), col="red")
  abline(v=-log(FoldChangeThreshold,2), col="red")
  abline(h=log(SignificanceThreshold,10)*-1, col="green")
  p_names <- row.names(ProtQuantMatrix_rn)
  proteinNamesPValFCs <- cbind(as.character(p_names), pValueVector, log2FCvector)
  return(proteinNamesPValFCs)
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

#CHECKED
# VirginzeMatrix
VirginzeMatrixWithNAs = function(ProtQuantMatrix_rn){
  NumSamples <- ncol(ProtQuantMatrix_rn)
  QMonly <- ProtQuantMatrix_rn
  #replace zeros with NA
  QMonly[QMonly==0] <- NA
  ProtQuantMatrix_rn_virg<-QMonly
  colnames(ProtQuantMatrix_rn_imputed) = colnames(ProtQuantMatrix_rn)  
  return(ProtQuantMatrix_rn_virg)
}


#Intensity_grp1vs2
PlotMeanGrp1vs2 = function(ProtQuantMatrix_rn, SignificanceThreshold){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
  }
  
  sigThreshold <- SignificanceThreshold
  bool_pVector<- vector()
  bool_pVector <- pValueVector<sigThreshold
  
  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  
  #xyHplot <- plot(log(grp1_means,10)~log(grp2_means,10), col=as.factor(bool_pVector))
  xyHplot <- plot(grp1_means~grp2_means, col=as.factor(bool_pVector))
  xyHplot <- abline(0,1, col="grey")
}



#Intensity_grp1vs2
PlotMeanGrp1vs2LogScale = function(ProtQuantMatrix_rn, SignificanceThreshold){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
  }
  
  sigThreshold <- SignificanceThreshold
  bool_pVector<- vector()
  bool_pVector <- pValueVector<sigThreshold
  
  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  
  #xyHplot <- plot(log(grp1_means,10)~log(grp2_means,10), col=as.factor(bool_pVector))
  xyHplot <- plot(grp1_means~grp2_means, col=as.factor(bool_pVector), log="xy")
  xyVolcanoPlot <- abline(0,1, col="grey")
}

#SmoothIntensity_grp1vs2
SmoothPlotMeanGrp1vs2LogScale = function(ProtQuantMatrix_rn, SignificanceThreshold){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
  }
  
  sigThreshold <- SignificanceThreshold
  bool_pVector<- vector()
  bool_pVector <- pValueVector<sigThreshold
  
  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  
  #xyHplot <- plot(log(grp1_means,10)~log(grp2_means,10), col=as.factor(bool_pVector))
  xyHplot <- smoothScatter(grp1_means~grp2_means, col=as.factor(bool_pVector))
  xyHplot <- abline(0,1, col="grey")
}


#GroupMeans
ReturnGroupMeans = function(ProtQuantMatrix_rn){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustTtest.returnpValue(grp_1[i,],grp_2[i,])
  }
  
  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  outVecBounded <- vector()
  outVecBounded <- cbind(row.names(ProtQuantMatrix_rn), grp1_means, grp2_means, pValueVector)  
  return(outVecBounded)
}

#GroupMeansOptionallyPaired
ReturnGroupMeansOptionallyPairedTest = function(ProtQuantMatrix_rn, Bool_pairedOrNot=FALSE){
  reps <- ncol(ProtQuantMatrix_rn)/2
  grp_1 <- ProtQuantMatrix_rn[,1:reps]
  grp_2 <- ProtQuantMatrix_rn[,(reps+1):ncol(ProtQuantMatrix_rn)]
  #ttest
  pValueVector <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    pValueVector[i] <- RobustPairedOrUnpairedTtest.returnpValue(grp_1[i,],grp_2[i,], pairedTRUEorFalse = Bool_pairedOrNot)
  }
  
  #mean vectors
  grp1_means <- vector()
  grp2_means <- vector()
  for (i in 1:nrow(ProtQuantMatrix_rn)) {
    grp1_means[i] <- rowMeans(grp_1[i,], na.rm=TRUE)
    grp2_means[i] <- rowMeans(grp_2[i,], na.rm=TRUE)
  }
  outVecBounded <- vector()
  outVecBounded <- cbind(row.names(ProtQuantMatrix_rn), grp1_means, grp2_means, pValueVector)  
  return(outVecBounded)
}



#CHECKED
#HeatMap on Proteins
DoHeatMapForProteins = function(ProtQuantMatrix_rn, ...){
  PQmat <- na.omit(as.matrix(ProtQuantMatrix_rn))
  heatmap.2(PQmat,margin=c(10,10),trace="none", ...)
}



#Normalization
# CHECKED
NormalizeTotalSumOfPQMatrix = function(ProtQuantMatrix_rn){
  maxSum <- max(apply(na.omit(ProtQuantMatrix_rn[,1:ncol(ProtQuantMatrix_rn)]),2,sum))
  ScaleFactorssum <- 1/(apply(na.omit(ProtQuantMatrix_rn[,1:ncol(ProtQuantMatrix_rn)]),2,sum)/maxSum)
  write.table(ScaleFactorssum, "appliedScaleFactors_SumNormalization.txt",row.names=TRUE,col.names=FALSE)
  nPQmatrix <- ProtQuantMatrix_rn
  for (i in 1:ncol(ProtQuantMatrix_rn)) {
    nPQmatrix[,i]<-ProtQuantMatrix_rn[,i]*ScaleFactorssum[i]
  }
  return(nPQmatrix)
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

#Actin normalization
NormalizePQMatrixWithHouseKeepingProtein = function(ProtQuantMatrix_rn, HousekeeperProtein){
  HouseKeeperLine <- agrep(HousekeeperProtein, row.names(ProtQuantMatrix_rn))
  MaxValueOfHK <- max(ProtQuantMatrix_rn[HouseKeeperLine,])
  ScaleFactorsHK <- 1/(ProtQuantMatrix_rn[HouseKeeperLine,]/MaxValueOfHK)
  write.table(ScaleFactorsHK, "appliedScaleFactors_HouseKeeperProt.txt",row.names=TRUE, col.names=FALSE)
  nPQmatrix <- ProtQuantMatrix_rn
  for (i in 1:ncol(ProtQuantMatrix_rn)) {
    nPQmatrix[,i]<-ProtQuantMatrix_rn[,i]*as.numeric(ScaleFactorsHK[i])
  }
  return(nPQmatrix)
}

#ProtQuantMatrix_rn <- n_i_dat
#SignificanceThreshold <- 0.05
#boolTrafoHere = TRUE
#FCcutoff <- 1.5
MAplotWith2grpsColorSignificants = function(ProtQuantMatrix_rn, SignificanceThreshold, FCcutoff, boolTrafoHere){
  pValVec <- vector()  
  TwogrpTest <- vector()  
  TwogrpTest<- Do2grpTtestOnMatrixAndBHcorrReturnAllInternalTrafo(ProtQuantMatrix_rn, bool_TrafoHere=boolTrafoHere)
  pValVec <- TwogrpTest[,2]
  isSignific <- vector()
  isSignific <- pValVec<SignificanceThreshold
  WhereIsNA <- vector()
  WhereIsNA <- is.finite(DoLog2FCCalc(ProtQuantMatrix_rn))
  logRatio <- vector()
  logRatio <- DoLog2FCCalc(ProtQuantMatrix_rn)
  all_means <- vector()
  all_means <- rowMeans(matrix(as.numeric(ReturnGroupMeans(ProtQuantMatrix_rn)[,2:3]),ncol=2))
  plot(logRatio[WhereIsNA] ~ sqrt(all_means[WhereIsNA]), ylim=c(-10,10), main="Bland-Altman(MA) plot", xlab="sqrt of mean abundance", ylab="log2Ratio",col=(as.numeric(isSignific)+1), pch=".")
  points(logRatio[WhereIsNA] ~ sqrt(all_means[WhereIsNA]), cex=0.2, col=(as.numeric(isSignific)+1))
  abline(h=median(logRatio), col="lightblue", lty=2, lwd=2)
  abline(h=0, col="grey")
  abline(h=log(FCcutoff,2), col="blue")
  abline(h=-log(FCcutoff,2), col="blue")
}

#ProtQuantMatrix_rn <- dat
ReturnValuesForMAplotWith2grps = function(ProtQuantMatrix_rn){
  pValVec <- vector()  
  pValVec<-Do2grpTtestOnMatrix(ProtQuantMatrix_rn)
  WhereIsNA <- vector()
  WhereIsNA <- is.finite(DoLog2FCCalc(ProtQuantMatrix_rn))
  logRatio <- vector()
  logRatio <- DoLog2FCCalc(ProtQuantMatrix_rn)
  all_means <- vector()
  all_means <- rowMeans(matrix(as.numeric(ReturnGroupMeans(ProtQuantMatrix_rn)[,2:3]),ncol=2))
  outVectorWithAllin <- vector()
  outVectorWithAllin <- cbind(as.character(ProtQuantMatrix_rn[WhereIsNA,1]), logRatio[WhereIsNA], all_means[WhereIsNA], pValVec[WhereIsNA])
  return(outVectorWithAllin)
}


#RawProtQuantMatrixRow <- quantDat[2577,]
#NormImputedProtQuantMatrixRow <- n_i_dat[2577,]
#Do it also for only 2 groups since it makes a lot of sense to see where the imputed values are..
DoBoxplotWithBeeswarmPointsAndRecalcAOVIndicateImputesFor2GroupsProperOrganized = function(NormImputedProtQuantMatrixRow, RawProtQuantMatrixRow, GroupingFactors){
  VectorOfImputes <- RawProtQuantMatrixRow!=0
  protValues <- as.double(NormImputedProtQuantMatrixRow)
  aovFit <- aov(asinh(protValues) ~ GroupingFactors)
  s <- summary(aovFit)
  boxplot(asinh(protValues)~ GroupingFactors, outline=F, main=rownames(NormImputedProtQuantMatrixRow), col="#ededed", xlab="Groups", ylab="arcsinh(ProteinAbundance)",sub = paste("t-test (pValue):", round(s[[1]][["Pr(>F)"]][1], 6)))
  colVec <- (1/(as.integer(VectorOfImputes)+1)*2)
  #swarm it
  beeswarm(asinh(protValues)~ GroupingFactors, pwcol=colVec, pwpch=colVec, lwd=15,add=TRUE)
  legend("bottomleft", legend = c("measured values", "imputed values"), pch=1:2, col=1:2)
}


#unchecked .. 
#x[!is.na(x)]
RemoveZerosFromProtMatrixNKeepRow = function(ProtQuantMatrix_rn){
  n <- ncol(ProtQuantMatrix_rn)
  columns <- c(1:n)
  ProtQuantMatrix_rn[ProtQuantMatrix_rn[,columns]==0] <- NA
  return(ProtQuantMatrix_rn)
} 


############################
#      Visualize Multigroup 
############################

#NormImputedProtQuantMatrixRow <- TopProteinLine
#RawProtQuantMatrixRow <- TopProteinRAWline
#GroupingFactors <- fGroups
#rm(NormImputedProtQuantMatrixRow, RawProtQuantMatrixRow, GroupingFactors)


DoBoxplotWithBeeswarmPointsAndRecalcAOVIndicateImputes = function(NormImputedProtQuantMatrixRow, RawProtQuantMatrixRow, GroupingFactors){
  VectorOfImputes <- RawProtQuantMatrixRow!=0
  protValues <- as.double(NormImputedProtQuantMatrixRow)
  aovFit <- aov(asinh(protValues) ~ GroupingFactors)
  s <- summary(aovFit)
  #mar <- c(1,1,1,1)
  #plot it   "#00000050"
#  boxplot(asinh(protValues)~ GroupingFactors, outline=F, main=rownames(NormImputedProtQuantMatrixRow), col=rainbow(length(unique(GroupingFactors))), xlab="Groups", ylab="arcsinh(ProteinAbundance)",sub = paste("ANOVA Pr(>F):", round(s[[1]][["Pr(>F)"]][1], 4)))
  boxplot(asinh(protValues)~ GroupingFactors, outline=F, main=rownames(NormImputedProtQuantMatrixRow), col="#ededed", xlab="Groups", ylab="arcsinh(ProteinAbundance)",sub = paste("ANOVA (pValue):", round(s[[1]][["Pr(>F)"]][1], 6)))
  colVec <- (1/(as.integer(VectorOfImputes)+1)*2)
  #swarm it
  beeswarm(asinh(protValues)~ GroupingFactors, pwcol=colVec, pwpch=colVec, lwd=15,add=TRUE)
  legend("bottomleft", legend = c("measured values", "imputed values"), pch=1:2, col=1:2)
}






