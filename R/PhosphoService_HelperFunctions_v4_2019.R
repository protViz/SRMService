# 2018-09-04: First running version
# 2019-01-09: Refactoring and better coding style to get it accepted in SRMService ;)
# 2019-09-03: Commiting PhosphoService functions into github

#' prepare smaller data.frame from original MaxQuant PhosphoSTYSite file and reduce to columns that are useful and filter with Localizations probability
#' @param phosphoSTYSitetableFromMQ original MaxQuant Phospho_STY_file
#' @param LocalizationProbabilityFilter Threshold on Localization.prob
#' @param min_nonNAinTotal minimum number of valid values (non NAs)
#' @return data.frame with less columns and rows
#' @export
#' @importFrom stringr str_count
#'
prepareSmallerMatrixFromPhospohSTY <- function(phosphoSTYtableFromMQ, LocalizationProbabilityFilter = 0.75,
                                               min_nonNAinTotal = 3) {
  pDat <- read.table(phosphoSTYtableFromMQ, header = TRUE, sep="\t", quote="")
  IntensityColumnsWithMultiplicity_idx <- grep(pattern = "___", x = colnames(pDat))
  qMat <- pDat[, IntensityColumnsWithMultiplicity_idx]
  # Cut PhosphoSTY down to useful columns
  includeList <- c("Localization.prob", "Protein", "Leading.proteins", "Positions.within.proteins","Protein.group.IDs",
                   "PEP", "Score", "Delta.score", "Position.in.peptide", "Phospho..STY..Probabilities", "Charge",
                   "Number.of.Phospho..STY.","Amino.acid", "Peptide.window.coverage","Sequence.window",
                   "Modification.window","Reverse", "Potential.contaminant")
  NonQMat <- pDat[,includeList]
  # Add column w naked sequence
  NonQMat$nakedSequence <- .extractNakedPepSequenceUsing2ColumnsFromMQ(seqWindowColumn = NonQMat$Sequence.window,
                                                                      PepWindowColumn = NonQMat$Peptide.window.coverage)
  # Go log2
  lg2Mat <- log2(qMat)
  # filter for localization probability
  bool_locP <- NonQMat$Localization.prob > LocalizationProbabilityFilter
  # before going: Expand Site Table
  perseusLikeMat <- cbind(NonQMat[bool_locP,], lg2Mat[bool_locP,])
  # Expand Site table
  gatheredMat <- .expandSiteTableFunction_forMQSiteTable(startMat = perseusLikeMat)

  # filter for too many NAs
  intMat <- gatheredMat[,grep(pattern = "Intensity", x = colnames(gatheredMat))]
  intMat[!is.finite(as.matrix(intMat))] <- NA

  # some filtering where we more than 32 - 3 NAs == more than 29!
  bool_GoodToKeep <- rowSums(!is.na(intMat[,])) > min_nonNAinTotal

  # here we generate the relevant matrix
  intMat_good <- intMat[bool_GoodToKeep,]
  nonQ_good <- gatheredMat[bool_GoodToKeep,-c(grep(pattern = "Intensity", x = colnames(gatheredMat)))]
  goodMatrixRearranged <- cbind(nonQ_good, intMat_good)

  # Prepare the Matrix for the merge and quantitative stuff
  myProtIdx <- 1:nrow(goodMatrixRearranged)
  myProteinNames <- paste(myProtIdx, "-", goodMatrixRearranged$Protein, "-",goodMatrixRearranged$PepModSeqProbabilities,
                          "-", goodMatrixRearranged$multiplicity, sep="")
  goodMatrixRearranged$myProteinNames <- myProteinNames
  return(goodMatrixRearranged)
}




# this undocumented function is used to blow up the matrix for the multiplicity in PhosphoSTY site table
.expandSiteTableFunction_forMQSiteTable = function(startMat){
  int_names <- colnames(startMat)[grep(pattern = "Intensity", x = colnames(startMat))]
  RawFilenameVector <- vector()
  for (i in 1:length(int_names)) {
    RawFilenameVector[i] <- stringr::str_split(string = int_names[i], pattern = "___")[[1]][1]
  }

  # the first element is "Intensity" from all rawFiles -> start from 2:x
  #length(unique(RawFilenameVector))
  allRFNuni <- unique(RawFilenameVector)[2:length(unique(RawFilenameVector))]
  int_idx <- grep(pattern = "Intensity", x = colnames(startMat))
  infoMat <- startMat[,-int_idx]
  outputMat <- matrix()
  for (i in 1:length(allRFNuni)) {
    idx_oi <- grep(x = colnames(startMat), pattern = allRFNuni[i])
    #matForReshaping <- cbind(infoMat, startMat[,idx_oi])
    #
    tblMat <- tbl_df(cbind(infoMat, startMat[,idx_oi]))
    #colnames(tblMat)
    # we have to add some more attributes that might be important later
    colnames(tblMat) <- c("LocProb","Protein", "LeadingProteins", "PositionsINproteins","ProteinGroupIDs", "PEP",
                          "Score", "DeltaScore", "PositionInpeptide", "PepModSeqProbabilities", "z","numPhos", "modAA",
                          "PepWindowCoverage", "SeqWindow", "modWindow",  "Rev", "Cont","nakedSequence","single",
                          "double", "triple")
    partForGathered <- tblMat %>% gather(multiplicity, RawFileIntensityOfInterest, -c(LocProb, Protein,LeadingProteins,
                                                                                     PositionsINproteins,
                                                                                     ProteinGroupIDs, PEP, Score,
                                                                                     DeltaScore, PositionInpeptide,
                                                                                     PepModSeqProbabilities, z, numPhos,
                                                                                     modAA, PepWindowCoverage,
                                                                                     SeqWindow, modWindow, Rev, Cont,
                                                                                     nakedSequence))
    # rename cols
    myCn <- colnames(tblMat)[1:(length(colnames(tblMat))-3)]
    colnames(partForGathered) <- c(myCn ,"multiplicity",allRFNuni[i])
    outputMat <- cbind(outputMat, partForGathered)
  }
  # get rid of repeated columns
  int_idx2 <- grep(pattern = "Intensity", x = colnames(outputMat))
  endMeta_idx <- int_idx2[1]-1
  #
  outputMatRedefined <- outputMat[, c(2:endMeta_idx, int_idx2)]
  return(outputMatRedefined)
}

# this undocumented function is used to extract the naked peptide sequence from the pepwindow column
.extractNakedPepSequenceUsing2ColumnsFromMQ <- function(seqWindowColumn, PepWindowColumn) {
  myWindow <- seqWindowColumn
  myWindowCov <- PepWindowColumn
  myPepLengths <- str_count(myWindowCov, "P")
  myStarts <- vector(length=length(myWindowCov))
  for (i in 1:length(myWindowCov)) {
    myStarts[i] <- regexpr(pattern = "P", text = myWindowCov[i], useBytes = FALSE)[1]
  }
  myStops <- myStarts + myPepLengths
  myNakedPeps <- substr(x = myWindow,start = myStarts, stop = myStops)
  return(myNakedPeps)
}

#' This function builds up the Annotations from the raw-file names if they follow the naming conventions
#' @param proteinMAT proteinMatrix object as read in
#' @param numericSeperatorForFilenames id_x where the conditions are encoded in the file names
#' @return annotation table
#' @export
#'
buildUpAnnotationFileFromColumnNames <- function(proteinMAT, numericSeperatorForFilenames = 5){
  # all raw files in protein groups (select here for proper 2grp)
  rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(proteinMAT),value=T) )
  # number 4 splits in treated vs untreated
  condition <- quantable::split2table(rawF)[,numericSeperatorForFilenames]
  #
  annotationTab <- data.frame(Raw.file = rawF,
                              Condition = condition,
                              BioReplicate = paste("X",1:length(condition),sep=""),
                              Run = 1:length(condition),
                              IsotopeLabelType = rep("L",length(condition)),
                              stringsAsFactors = F)
  return(annotationTab)
}


#' This function returns the result of the 2grp analysis as tables as well as object
#' @param proteinMat quantitative matrix (can alos be peptide matrix)
#' @param pathNFileForAnnotation working directory where the annotation file is found
#' @param mxNAs maximum number of NAs allowed before rows are filtered-out
#' @param nrPep minimum number of peptides
#' @param DenomReference reference that is used as denominator in cases where ratios (or logratios) are used (e.g. volcano plot)
#' @param expName name for the experiment
#' @param resultDirectory directory generated on the fly where the annotation table is saved
#' @param lg2FCthreshold threshold for the absolute log2 fold-change (symmetric in both directions)
#' @param qvThreshold threshold for qValue (adjust.p.value)
#' @param Pairing boolean TRUE/FALSE
#' @return list with results of a twoGroup analysis
#' @export
#'
get2GroupAnalysisFromExistingAnnotationFileNoRender <- function(proteinMat, pathNFileForAnnotation,
                                                                             mxNAs, nrPep, DenomReference,
                                                                             expName = "GlobalProtein_comp1",
                                                                             resultDirectory = "myResultDirectory",
                                                                             lg2FCthreshold = qfoldchange,
                                                                             qvThreshold = qvalueThreshold,
                                                                             Pairing=FALSE){
  Experimentname = expName
  resultdir <- resultDirectory
  dir.create(resultdir)
  # read in annotation which are fixed
  annotationInFn <- read.table(file = pathNFileForAnnotation, header=TRUE)
  # Anno file name
  annoFN <- paste(file.path(resultdir, paste("AnnotationFile_fn_", Experimentname,".txt", sep="")))
  write.table(annotationInFn, file=annoFN)
  # make it more robust.. copy  particular Rnw files directly from package here!
  reference=DenomReference
  #print("A")
  grp2 <- Grp2Analysis(annotationInFn, Experimentname, maxNA=mxNAs, nrPeptides=nrPep, reference=reference)
  #print("B")
  grp2$setMQProteinGroups(proteinMat)
  grp2$setQValueThresholds(qvalue = qvThreshold, qfoldchange = lg2FCthreshold)
  grp2$annotation_$Condition <- factor(grp2$annotation_$Condition,
                                       levels = c(grp2$reference, setdiff(grp2$conditions,grp2$reference )))

  # Here set pairing!!
  IsPaired <- Pairing

  if(IsPaired){
    x1 <- model.matrix(~ Condition + BioReplicate, grp2$annotation_)
    grp2$setModelMatrix(x1)
    myResultTable <- grp2$getResultTable()
  }else{
    x1 <- model.matrix(~ Condition, grp2$annotation_)
    grp2$setModelMatrix(x1)
    myResultTable <- grp2$getResultTable()
  }
  return(list(myResultTable = myResultTable, grp2 = grp2))
}



# Do the NtoCplot
# this undocumented function is used to plot the phospho peptides from NtoC-term with indicated log2-fold-change-bars
# the function is used inside another documented function
.plotProteinNPhosphoPeptidesNtoC <- function(protName, globalProtFC, POI_tableProtein, POI_tablePeptide,
                                            protColor = "blue", protWtdh = 4, maxX = 1200, protNCoffset = 40,
                                            pepSigStarOffset = 0.1, pepSigStarSize = 2, qModThreshold = 0.3,
                                            pModThreshold = 0.1, greyLineLog2Threshold = absLg2FCThreshold) {
  myYaxisLimiter <- max(ceiling(max(abs(POI_tablePeptide$log2FC))), ceiling(max(abs(POI_tableProtein$log2FC))),
                        ceiling(max(abs(POI_tablePeptide$pseudoLog2FC))), na.rm = TRUE)
  miny <- -myYaxisLimiter
  maxy <- myYaxisLimiter
  # get only one (the first if multiple) out
  Phos_PositionsInProteins <- vector(length=nrow(POI_tablePeptide))
  for (i in 1:nrow(POI_tablePeptide)) {
    if(str_count(string = POI_tablePeptide$PositionsINproteins[i] , pattern = ";") == 0) {
      Phos_PositionsInProteins[i] <- POI_tablePeptide$PositionsINproteins[i]
    } else {
      Phos_PositionsInProteins[i] <- as.numeric(stringr::str_split(string = POI_tablePeptide$PositionsINproteins[i],
                                                          pattern = ";")[[1]][1])
    }
  }
  scaleX <- 1/max(Phos_PositionsInProteins)*1200
  POI_tablePeptide$PositionsINproteins <- Phos_PositionsInProteins
  phosphoPlotTitle <- paste("Prot: ",protName,"\n Peptides from: ", POI_tablePeptide$Protein[1],"# phospho = ",
                            length(POI_tablePeptide$log2FC), sep=" ")
  # Do the plot
  plot(c(-100,0,1400, 0,0), c(globalProtFC,0,0,miny, maxy), pch=".", col="white", main=phosphoPlotTitle, ylab="log2FC",
       xlab="Full length protein", xaxt='n')
  axis(side = 1,at = c(-100,0,1200),labels = c("protFC","0","100%"))
  #grey box for prot FC
  rect(xleft = -150,ybottom = -100,xright = -50,ytop = 100, col = "lightgrey",border = FALSE)
  if (globalProtFC == 0) {
    points(-100,0,pch="0", col="orange", cex=3)
  } else {
    segments(-100,0,-100, globalProtFC, col=protColor, lwd=protWtdh)
  }
  #segments(-120,0,-80,0,col="black", lwd=protWtdh)
  segments(0,0,maxX,0, col="black", lwd=1.5)
  points(-protNCoffset,0,pch="N")
  points(maxX+protNCoffset,0, pch="C")
  # plot all peptides
  for (i in 1:length(POI_tablePeptide$PositionsINproteins)) {
    # check if one sider !!
    if (is.na(POI_tablePeptide$log2FC[i])) {
      segments(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX, 0,
               as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX,
               POI_tablePeptide$pseudoLog2FC[i], lwd=1, lty=2, col="orange")

      if(POI_tablePeptide$pseudoLog2FC[i] > 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                      scaleX, POI_tablePeptide$pseudoLog2FC[i]+pepSigStarOffset,
                                                      pch="x", cex=pepSigStarSize,
                                                      col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

      if(POI_tablePeptide$pseudoLog2FC[i] < 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                      scaleX, POI_tablePeptide$pseudoLog2FC[i]-pepSigStarOffset,
                                                      pch="x", cex=pepSigStarSize,
                                                      col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

    } else {
      # no one sider peptide
      segments(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX, 0,
               as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX, POI_tablePeptide$log2FC[i],
               col=as.numeric(as.factor(POI_tablePeptide$modAA[i])), lwd=1)

      if (POI_tablePeptide$q.mod[i] <= qModThreshold) {
        if(POI_tablePeptide$log2FC[i] > 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]+pepSigStarOffset, pch="*",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

        if(POI_tablePeptide$log2FC[i] < 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]-pepSigStarOffset, pch="*",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

      } else if (POI_tablePeptide$p.mod[i] < pModThreshold) {
        if(POI_tablePeptide$log2FC[i] > 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]+pepSigStarOffset, pch="+",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

        if(POI_tablePeptide$log2FC[i] < 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]-pepSigStarOffset, pch="+",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))
      }
    }
  }
  legend("bottomright", legend = unique(POI_tablePeptide$modAA),
         col=unique(as.numeric(as.factor(POI_tablePeptide$modAA))),
         text.col = unique(as.numeric(as.factor(POI_tablePeptide$modAA))), lwd=1)
  legend("topright", legend = c("PseudoFoldChange", "Significant q-Value", "Significant p-Value"),
         col = c("Orange", "black", "black"), pch=c("x","*", "+"), cex = 0.7)
  abline(h=greyLineLog2Threshold, col="grey", lty=2)
  abline(h=-greyLineLog2Threshold, col="grey", lty=2)
}


# Do the NtoCplot
# this undocumented function is used to plot the phospho peptides from NtoC-term with indicated log2-fold-change-bars
# the function is used inside another documented function -> the 2019 function should be used from 2019 onwards
.plotProteinNPhosphoPeptidesNtoC_2019 <- function(protName, globalProtFC, POI_tableProtein, POI_tablePeptide,
                                                 protColor = "blue", protWtdh = 4, maxX = 1200, protNCoffset = 40,
                                                 pepSigStarOffset = 0.1, pepSigStarSize = 2, qModThreshold = 0.3,
                                                 pModThreshold = 0.1, greyLineLog2Threshold = absLg2FCThreshold) {
  myYaxisLimiter <- max(ceiling(max(abs(POI_tablePeptide$log2FC))), ceiling(max(abs(POI_tableProtein$log2FC))),
                        ceiling(max(abs(POI_tablePeptide$pseudoLog2FC))), na.rm = TRUE)
  miny <- -myYaxisLimiter
  maxy <- myYaxisLimiter
  # get only one (the first if multiple) out
  Phos_PositionsInProteins <- vector(length=nrow(POI_tablePeptide))
  for (i in 1:nrow(POI_tablePeptide)) {
    if(str_count(string = POI_tablePeptide$PositionsINproteins[i] , pattern = ";") == 0) {
      Phos_PositionsInProteins[i] <- POI_tablePeptide$PositionsINproteins[i]
    } else {
      Phos_PositionsInProteins[i] <- as.numeric(stringr::str_split(string = POI_tablePeptide$PositionsINproteins[i],
                                                          pattern = ";")[[1]][1])
    }
  }
  scaleX <- 1/max(Phos_PositionsInProteins)*1200
  POI_tablePeptide$PositionsINproteins <- Phos_PositionsInProteins
  phosphoPlotTitle <- paste("Prot: ",protName,"\n Peptides from: ", POI_tablePeptide$Protein[1],"# phospho = ",
                            length(POI_tablePeptide$log2FC), sep=" ")
  # Do the plot
  plot(c(-100,0,1400, 0,0), c(globalProtFC,0,0,miny, maxy), pch=".", col="white", main=phosphoPlotTitle, ylab="log2FC",
       xlab="Full length protein", xaxt='n')
  axis(side = 1,at = c(-100,0,1200),labels = c("protFC","0","100%"))
  #grey box for prot FC
  rect(xleft = -150,ybottom = -100,xright = -50,ytop = 100, col = "lightgrey",border = FALSE)
  if (globalProtFC == 0) {
    points(-100,0,pch="0", col="orange", cex=3)
  } else {
    segments(-100,0,-100, globalProtFC, col=protColor, lwd=protWtdh)
  }
  #segments(-120,0,-80,0,col="black", lwd=protWtdh)
  segments(0,0,maxX,0, col="black", lwd=1.5)
  points(-protNCoffset,0,pch="N")
  points(maxX+protNCoffset,0, pch="C")
  # plot all peptides
  for (i in 1:length(POI_tablePeptide$PositionsINproteins)) {
    # check if one sider !!
    if (is.na(POI_tablePeptide$log2FC[i])) {
      segments(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX, 0,
               as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX,
               POI_tablePeptide$pseudoLog2FC[i], lwd=1, lty=2, col="orange")

      if(POI_tablePeptide$pseudoLog2FC[i] > 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                        scaleX, POI_tablePeptide$pseudoLog2FC[i]+pepSigStarOffset,
                                                      pch="x", cex=pepSigStarSize,
                                                      col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

      if(POI_tablePeptide$pseudoLog2FC[i] < 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                        scaleX, POI_tablePeptide$pseudoLog2FC[i]-pepSigStarOffset,
                                                      pch="x", cex=pepSigStarSize,
                                                      col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

    } else {
      # no one sider peptide
      segments(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX, 0,
               as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i]))*scaleX, POI_tablePeptide$log2FC[i],
               col=as.numeric(as.factor(POI_tablePeptide$modAA[i])), lwd=1)

      if (POI_tablePeptide$adj.P.Val[i] <= qModThreshold) {
        if(POI_tablePeptide$log2FC[i] > 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]+pepSigStarOffset, pch="*",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

        if(POI_tablePeptide$log2FC[i] < 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]-pepSigStarOffset, pch="*",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

      } else if (POI_tablePeptide$P.Value[i] < pModThreshold) {
        if(POI_tablePeptide$log2FC[i] > 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]+pepSigStarOffset, pch="+",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))

        if(POI_tablePeptide$log2FC[i] < 0) points(as.numeric(as.character(POI_tablePeptide$PositionsINproteins[i])) *
                                                    scaleX, POI_tablePeptide$log2FC[i]-pepSigStarOffset, pch="+",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(as.factor(POI_tablePeptide$modAA[i])))
      }
    }
  }
  legend("bottomright", legend = unique(POI_tablePeptide$modAA),
         col=unique(as.numeric(as.factor(POI_tablePeptide$modAA))),
         text.col = unique(as.numeric(as.factor(POI_tablePeptide$modAA))), lwd=1)
  legend("topright", legend = c("PseudoFoldChange", "Significant q-Value", "Significant p-Value"),
         col = c("Orange", "black", "black"), pch=c("x","*", "+"), cex = 0.7)
  abline(h=greyLineLog2Threshold, col="grey", lty=2)
  abline(h=-greyLineLog2Threshold, col="grey", lty=2)
}


#' prepare generic intensityMatrix from the subsetted and handled matrix before
#' @param myGoodMatrix the quant matrix prepared with the prepareSmallerMatrixFromPhospohSTY function (note, intensities are log2 in the input and linear in the output)
#' @return data.frame with less columns and rows
#' @export
#' @examples
#' prepareGenericIntensityMatrixFromGoodMatrixForTesting(SRMService::myMat_good)
#'
prepareGenericIntensityMatrixFromGoodMatrixForTesting <- function(myGoodMatrix) {
  fakeNrPeptides <- rep(3, nrow(myGoodMatrix))
  myProteinNames <- myGoodMatrix$myProteinNames
  fakeFastaHeader <- rep("No_Description", nrow(myGoodMatrix))
  relevant_intMat <- myGoodMatrix[,grep(pattern = "Intensity", x = colnames(myGoodMatrix))]
  relevant_intMatNonLog2 <- 2^relevant_intMat
  MatrixForQuantTesting <- data.frame(myProteinNames, relevant_intMatNonLog2,fakeNrPeptides, fakeFastaHeader)
  newHeaders <- c("Majority.protein.IDs", colnames(relevant_intMatNonLog2), "Peptides","Fasta.headers")
  # 2019-11-22: important fix to see Majority... as character and not as factors
  colnames(MatrixForQuantTesting) <- newHeaders
  MatrixForQuantTesting$Majority.protein.IDs <- as.character(MatrixForQuantTesting$Majority.protein.IDs)
  return(MatrixForQuantTesting)
}


#' merge statistics from 2grp with the additional information from the phosphoSTY file into one
#' @param myGoodMatrix the quant matrix prepared with the prepareSmallerMatrixFromPhospohSTY function
#' @param phosPepsQuantResults the table from the 2grp analysis with the statistics
#' @return merged data frame
#' @export
#' @examples
#' getDFwithMergedQuantResultsWithLargeMatrix(theGoodMatrix = SRMService::myMat_good, phosPepsQuantResults = SRMService::quantResPhosphoPep)
#'
getDFwithMergedQuantResultsWithLargeMatrix <- function(theGoodMatrix, phosPepsQuantResults) {
  # Prepare the Matrix for more
  myProtIdx <- 1:nrow(theGoodMatrix)
  myProteinNames <- paste(myProtIdx, "-", theGoodMatrix$Protein, "-",theGoodMatrix$PepModSeqProbabilities, "-",
                          theGoodMatrix$multiplicity, sep="")
  theGoodMatrix$myProteinNames <- myProteinNames
  FullResultDF <- merge(x = theGoodMatrix, y = phosPepsQuantResults, by.x = "myProteinNames", by.y = "ProteinName")
  return(FullResultDF)
}

#' Subset merged matrix based on thresholds for moderated p.Value and absolute log2 fold change
#' @param mergedResultMatrix result matrix with statistics from 2grp analysis and phosphoSTY information
#' @param my_pModThreshold threshold for moderated p.value
#' @param my_absLog2FCThreshold threshold absolute log2 fold-change
#' @return candidates (subset matching thresholds)
#' @export
#'
getCandidatesFromMergedMatrix <- function(mergedResultMatrix, my_pModThreshold, my_absLog2FCThreshold) {
  #bool_pModOK <- mergedResultMatrix$q.mod < my_qvThreshold
  bool_pModOK <- mergedResultMatrix$P.Value < my_pModThreshold
  bool_abslg2FCOK <- abs(mergedResultMatrix$log2FC) > my_absLog2FCThreshold
  # pseudo relevant (full NAs in one condition!)
  bool_sig_phosPeps <- vector(length=nrow(mergedResultMatrix))
  for (i in 1:nrow(mergedResultMatrix)) {
    # here we look for thresholds to be met -> the pseudo-significant ones are handled below
    try (if (bool_abslg2FCOK[i] == TRUE && bool_pModOK[i] == TRUE) bool_sig_phosPeps[i] <- TRUE, silent = TRUE)
  }
  #sum(bool_sig_phosPeps)
  # here our hot candidates
  onesideCandidates <- mergedResultMatrix[which(is.na(bool_pModOK)),]
  # here the significant candidates
  significantCandidates <-  mergedResultMatrix[which(bool_sig_phosPeps),]
  myCandidates <- rbind(significantCandidates, onesideCandidates)
  return(myCandidates)
}


#' Write txt files with unique protein names that can be used for ORA analysis
#' @param mergedResultMatrix result matrix with statistics from 2grp analysis and phosphoSTY information
#' @param my_pModThreshold threshold for moderated p.value
#' @param my_absLog2FCThreshold threshold absolute log2 fold-change
#' @param comparisonName Name for comparison used for the filenames of txt files
#' @param geneNamesInsteadOfAccessions boolean to inidicate parsing uniprot accession or gene names sp|accession|gene_name
#' @param isUniprotDB boolean to indicate if identifier is of type uniprot
#' @export
#' @importFrom stringr str_split
#'
generateTxtFilesForDownstreamAnalysisAndReturnCandidateMatrix_2019 <- function(comparisonName="TestComparison",
                                                                          mergedResultMatrix,
                                                                          my_pModThreshold = pModThreshold,
                                                                          my_absLog2FCThreshold = absLg2FCThreshold,
                                                                          geneNamesInsteadOfAccessions = FALSE,
                                                                          isUniprotDB = TRUE){
  bgFileName <- paste("Background_", comparisonName, ".txt",sep="")
  fgFileName <- paste("DifferentiallyRegulated_", comparisonName, ".txt",sep="")
  #fgFileName_up <- paste("UpRegulated_", comparisonName, ".txt",sep="")
  #fgFileName_down <- paste("DownRegulated_", comparisonName, ".txt",sep="")

  #BG
  if (isUniprotDB == TRUE) {
    all_accessions <- sort(unlist(lapply(stringr::str_split(string  = mergedResultMatrix$Protein, pattern  = "\\|"),FUN = "[",2)))
    all_geneNames <- sort(unlist(lapply(stringr::str_split(string  = mergedResultMatrix$Protein, pattern  = "\\|"),FUN = "[",3)))
  } else {
    all_accessions <- sort(unlist(lapply(stringr::str_split(string  = mergedResultMatrix$Protein, pattern  = "\\|"),FUN = "[",1)))
    all_geneNames <- sort(unlist(lapply(stringr::str_split(string  = mergedResultMatrix$Protein, pattern  = "\\|"),FUN = "[",1)))
  }
  # prepare for outputfiles
  unique_all_accessionsWithCounts <- rle(all_accessions)
  #length(unique_all_accessionsWithCounts$values)
  unique_all_geneNamesWithCounts <- rle(all_geneNames)
  #length(unique_all_geneNamesWithCounts$values)

  # Filter down
  # here variable names are changed! p.mod -> P.Value
#  bool_pModOK <- mergedResultMatrix$p.mod < my_pModThreshold
  bool_pModOK <- mergedResultMatrix$P.Value < my_pModThreshold
  bool_abslg2FCOK <- abs(mergedResultMatrix$log2FC) > my_absLog2FCThreshold
  # pseudo relevant (full NAs in one condition!)
  bool_sig_phosPeps <- vector(length=nrow(mergedResultMatrix))
  for (i in 1:nrow(mergedResultMatrix)) {
    # check if thresholds are met -> pseudo-significant ones are handled below
    try (if (bool_abslg2FCOK[i] == TRUE && bool_pModOK[i] == TRUE) bool_sig_phosPeps[i] <- TRUE, silent = TRUE)
  }

  # here our hot candidates
  onesideCandidates <- mergedResultMatrix[which(is.na(bool_pModOK)),]
  # here the significant candidates
  significantCandidates <-  mergedResultMatrix[which(bool_sig_phosPeps),]
  myCandidates <- rbind(significantCandidates, onesideCandidates)

  #FG
  if (isUniprotDB == TRUE) {
    sig_accessions <- rle(sort(unlist(lapply(stringr::str_split(string  = myCandidates$Protein, pattern  = "\\|"),FUN = "[",2))))
    sig_geneNames <- rle(sort(unlist(lapply(stringr::str_split(string  = myCandidates$Protein, pattern  = "\\|"),FUN = "[",3))))
  } else {
    sig_accessions <- rle(sort(unlist(lapply(stringr::str_split(string  = myCandidates$Protein, pattern  = "\\|"),FUN = "[",1))))
    sig_geneNames <- rle(sort(unlist(lapply(stringr::str_split(string  = myCandidates$Protein, pattern  = "\\|"),FUN = "[",1))))
  }


  if (geneNamesInsteadOfAccessions) {
    write.table(x = unique_all_geneNamesWithCounts$values, file = paste("GeneNames_", bgFileName,sep=""),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(x = sig_geneNames$values, file = paste("GeneNames_", fgFileName,sep=""), row.names = FALSE,
                col.names = FALSE, quote = FALSE)

    write.table(rbind(c("GeneName", "counts"),cbind(unique_all_geneNamesWithCounts$values,
                                                    unique_all_geneNamesWithCounts$lengths)),
                file = paste("GeneNames_wCounts_", bgFileName), sep="\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE)

    write.table(rbind(c("GeneName", "counts"),cbind(sig_accessions$values, sig_accessions$lengths)),
                file = paste("GeneNames_wCounts_", fgFileName), sep="\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE)

  } else {
    write.table(x = unique_all_accessionsWithCounts$values, file = paste("Acc_", bgFileName,sep=""), row.names = FALSE,
                col.names = FALSE, quote = FALSE)

    write.table(x = sig_accessions$values, file = paste("Acc_", fgFileName,sep=""), row.names = FALSE,
                col.names = FALSE, quote = FALSE)

    write.table(rbind(c("Acc", "counts"),cbind(unique_all_accessionsWithCounts$values,
                                               unique_all_accessionsWithCounts$lengths)),
                file = paste("Acc_wCounts_", bgFileName), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    write.table(rbind(c("Acc", "counts"),cbind(sig_accessions$values, sig_accessions$lengths)),
                file = paste("Acc_wCounts_", fgFileName), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}



#' Simply merge phosphoAnalysis with total proteinAnalysis
#' @param mergedResultMatrix result matrix with statistics from 2grp analysis and more information from PhosphoSTY
#' @param globalProteinQuantResults result matrix with statistics from 2grp analysis for total protein
#' @return merged result lists with doublications on protein side (y)
#' @export
#' @examples
#' combineAndReturnPhosphoMatrix(mergedResultMatrix = SRMService::res__merge_allPhosPep, globalProteinQuantResults = SRMService::quantResTotalProt)
#'
combineAndReturnPhosphoMatrix <- function(mergedResultMatrix, globalProteinQuantResults){
  combinedNmergedResultMatrix <- merge(x = mergedResultMatrix, y = globalProteinQuantResults, by.x = "Protein",
                                       by.y = "TopProteinName", all.x = TRUE)
  return(combinedNmergedResultMatrix)
}


# 2019-08-23: adapted for updated SRM service package
#' Generate PDFs with NtoCplots for all candidates (using all phospho peptides even if not significant)
#' @param globalNphosphoCombinedNResultMatrix combinded result matrix with phospho and total protein result
#' @param candidateMatrix contains the significant phospho peptides for which all proteins should be drawn
#' @param expName this name is used to label the pdf file that is generated with this function
#' @export
#'
generateNtoCProteinPDFsWithPhosphoPeptides_2019 <- function(globalNphosphoCombinedNResultMatrix, candidateMatrix, expName) {
  # some static values to fill NAs
  linValueFiller <- 1
  qModForOneSiderCandidates <- 0.001
  pdfFileName <- paste("SignificantProtein_NtoCplot_", expName, ".pdf", sep="")
  mySigProteinHits <- rle(as.vector(sort(candidateMatrix$Protein)))
  # Do open pdf here
  pdf(pdfFileName,10,10)
  for (j in 1:length(mySigProteinHits$values)) {
    rm(POI_matrix)
    POI <- mySigProteinHits$values[j]
    # Extract Relevant Lines from Full table again
    POI_matrix <- globalNphosphoCombinedNResultMatrix[which(globalNphosphoCombinedNResultMatrix$Protein == POI), ]
    if(is.na(POI_matrix$P.Value.y[1])) {
      message("Protein Globally NOT quantified")
      MyProteinName_poi <- "Protein_globally_NOT_quantified"
      POI_protFC <- 0
    } else {
      message("We got it also globally")
      POI_protFC <- mean(POI_matrix$log2FC.y)
      MyProteinName_poi <- unique(POI_matrix$Protein)
    }

    # also clean phosphoPeps if NAs in
    idx_split <- grep(x = colnames(POI_matrix), pattern = "nrPeptides.y")
    grpSize <- length(grep(x = colnames(POI_matrix)[1:idx_split], pattern = "raw"))/2
    grp1_idx <- grep(x = colnames(POI_matrix)[1:idx_split], pattern = "raw")[1:grpSize]
    grp2_idx <- grep(x = colnames(POI_matrix)[1:idx_split],
                     pattern = "raw")[(grpSize+1):length(grep(x = colnames(POI_matrix)[1:idx_split], pattern = "raw"))]

    # Fill NAs with numbers
    # lin value filler
    for (f in 1:nrow(POI_matrix)) {
      for (ff in 1:length(grep(x = colnames(POI_matrix)[1:idx_split], pattern = "raw"))) {
        if (is.na(POI_matrix[f,grep(x = colnames(POI_matrix)[1:idx_split], pattern = "raw")[ff]])) {
          POI_matrix[f,grep(x = colnames(POI_matrix)[1:idx_split], pattern = "raw")[ff]] <- linValueFiller
        }
      }
    }

    # Do the plot
    # split matrix in protein part and peptide part
    idx_split <- grep(x = colnames(POI_matrix), pattern = "nrPeptides.y")
    POI_phosPeps<- POI_matrix[,1:(idx_split-1)]
    POI_globProts<- POI_matrix[,(idx_split-1):ncol(POI_matrix)]

    # PseudoFC are usually around 20!! -> Lets divide by two so to have a bit a better overview
    POI_phosPeps$pseudoLog2FC <-(log2(rowMeans(POI_phosPeps[,grp1_idx]) + 1) -
                                   log2(rowMeans(POI_phosPeps[,grp2_idx]) + 1))/2
    POI_phosPeps$adj.P.Val.x[is.na(POI_phosPeps$adj.P.Val.x)] <- qModForOneSiderCandidates
    # Do the plotting
    .plotProteinNPhosphoPeptidesNtoC_2019(protName = MyProteinName_poi, globalProtFC = POI_protFC,
                                    POI_tableProtein = POI_globProts, POI_tablePeptide = POI_phosPeps)
  }
  dev.off()
}




# Sequence logo plots
#' This functions generates a sequence logo from sequence windows (usually 31 chars long, PTMsite in the middle)
#' @param merged_phosphoPeptideResultMatrix combinded phospho result matrix with stats and info from PhosphoSTY
#' @param phosAcceptorAA Phospho acceptor site (in case of phospho: S|T|Y)
#' @param myLogoPlotTitle Title for the plot
#' @export
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom ggplot2 ggtitle
#'
doSequenceLogoPlotsFromCandidateMatrixForSinglePhosPeptides <- function(merged_phosphoPeptideResultMatrix,
                                                                        phosAcceptorAA = "S",
                                                                        myLogoPlotTitle = "Sequence Logo") {

  ProperLength <- median(nchar(as.vector(merged_phosphoPeptideResultMatrix$SeqWindow)))
  # count from how many we do the logo
  bool_ForPlotting <- merged_phosphoPeptideResultMatrix$multiplicity == "single" &
    merged_phosphoPeptideResultMatrix$modAA == phosAcceptorAA &
    nchar(as.vector(merged_phosphoPeptideResultMatrix$SeqWindow)) == ProperLength

  betterTitle <- paste(myLogoPlotTitle, " -> from: ",sum(bool_ForPlotting), " sequences", sep="")
  ggseqlogo(as.vector(merged_phosphoPeptideResultMatrix$SeqWindow[bool_ForPlotting]), seq_type='aa' ) +
    ggtitle(betterTitle)
}

