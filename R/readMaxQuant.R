.cast_maxquant_to_wide_glf <- function(d_long, aggregateFun=aggregateFun){
  data_w = dcast( Proteins + Modified.sequence + Charge ~ Raw.file, data=d_long, value.var='Intensity', fun.aggregate=aggregateFun, keep=TRUE)
  ## keep=TRUE : will keep the data.frame value as 1 even though there is no values for certain feature and certain run.
  ## when there is completely missing in certain feature and certain run, '1' will be filled. Therefore put NA instead of 1.
  data_w[data_w == 1] <- NA
  return(data_w)
}

.melt_maxquant_to_long_glf <- function(d_wide){
  data_l = melt(d_wide, id.vars=c('Proteins', 'Modified.sequence', 'Charge'))
  colnames(data_l)[colnames(data_l) %in% c("variable", "value")] <- c('Raw.file', 'Intensity')
  return(data_l)
}

.remove_feature_with_few <- function(x){
  count_measure = apply (x[, !(colnames(x) %in% c("Proteins", "Modified.sequence", "Charge"))], 1, function ( x ) length ( x[!is.na(x)] ) )
  remove_feature_name <- x[count_measure < 3, c("Proteins", "Modified.sequence", "Charge")]

  x$Feature <- paste(x$Proteins, x$Modified.sequence, x$Charge, sep="_")
  remove_feature_name$Feature <- paste(remove_feature_name$Proteins, remove_feature_name$Modified.sequence, remove_feature_name$Charge, sep="_")

  x <- x[-which(x$Feature %in% remove_feature_name$Feature), ]
  x <- x[, -ncol(x)]

  return(x)
}

################################################
### 1.1 remove contaminant, reverse proteinID
### Contaminant, Reverse column in evidence
.filterEvidence <- function(evidence){
  if(is.element("Contaminant", colnames(evidence)) & is.element("+", unique(evidence$Contaminant))){
    evidence <- evidence[-which(evidence$Contaminant %in% "+"), ]
  }

  if(is.element("Reverse", colnames(evidence)) & is.element("+", unique(evidence$Reverse))){
    evidence <- evidence[-which(evidence$Reverse %in% "+"),]
  }

  ### ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
  if(is.element("Only.identified.by.site", colnames(evidence))){
    evidence <- evidence[-which(evidence$Only.identified.by.site %in% "+"), ]
  }
  return(evidence)
}

.filterProteinGroups <- function(proteinGroups){
  ### first, remove contaminants
  if(is.element("Contaminant", colnames(proteinGroups)) & is.element("+",unique(proteinGroups$Contaminant))){
    proteinGroups <- proteinGroups[-which(proteinGroups$Contaminant %in% "+"), ]
  }

  if(is.element("Potential.contaminant", colnames(proteinGroups)) & is.element("+",unique(proteinGroups$Potential.contaminant))){
    proteinGroups <- proteinGroups[-which(proteinGroups$Potential.contaminant %in% "+"), ]
  }

  if(is.element("Reverse", colnames(proteinGroups)) & is.element("+",unique(proteinGroups$Reverse))){
    proteinGroups <- proteinGroups[-which(proteinGroups$Reverse %in% "+"), ]
  }

  ### ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
  if(is.element("Only.identified.by.site", colnames(proteinGroups))){
    proteinGroups <- proteinGroups[-which(proteinGroups$Only.identified.by.site %in% "+"), ]
  }
  return(proteinGroups)
}


#' MaxQtoMSstatsFormat
#' @param evidence MQ evidence.txt read using read.csv(sep="\t")
#' @param annotation  data.frame with columns Raw.file, Condition, BioReplicate, Run, IsotopeLabelType
#' @param proteinGroups provide MQ proteinGroups.txt
#' @param useUniquePeptide : (likely same as proteotypic) remove peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows : max or sum - when there are multiple measurements for certain feature and certain fun, use highest or sum of all.
#' @param fewMeasurements : if 1 or 2 measurements across runs per feature, 'remove' will remove those featuares. It can affected for unequal variance analysis.
#' @export
MQtoMSstatsFormat <- function(evidence,
                                annotation,
                                proteinGroups,
                                useUniquePeptide=TRUE,
                                summaryforMultipleRows=max,
                                fewMeasurements="remove",
                                removeMpeptides=TRUE){


  ### annotation.txt : Raw.file, Condition, BioReplicate, Run, (IsotopeLabelType)
  annot <- annotation


  evidence <- .filterEvidence( evidence )

  ################################################
  ### 1.1.2 matching proteinGroupID protein list

  ### need to check proteinGroupID in evidence and proteinGroup.txt the same
  ### 'id' in proteinGroups.txt vs 'Protein.group.IDs' in evidence
  ### possible to have some combination in Protein.group.IDs in evidence, such as 64;1274;1155;1273 instead of 64, 1274.. separately. combination of some ids seems not to be used for intensity
  ### 2015/02/03


  tempprotein <- .filterProteinGroups(proteinGroups)

  ### then take proteins which are included
  evidence <- evidence[which(evidence$Protein.group.IDs %in% unique(tempprotein$id)), ]

  ### then use 'protein.IDs' in proteinGroups.txt
  ### because if two 'proteins' in evidence.txt are used in one protein ID, need to use certain protein name in evidence.
  ### for example, protein.IDs in proteinGroups.txt are P05204;O00479. but, two 'proteins in evidence.txt, such as P05204;O00479, and P05204.

  tempname <- unique(tempprotein[,c("Protein.IDs", "id")])
  colnames(tempname) <- c("uniqueProteins", "Protein.group.IDs")

  evidence <- merge(evidence, tempname, by="Protein.group.IDs")

  evidence  <-  evidence[c("uniqueProteins", "Protein.group.IDs", "Sequence", "Modified.sequence", "Charge", "Raw.file", "Intensity", "Retention.time", "id")]

  colnames(evidence)[colnames(evidence) == "uniqueProteins"] <- "Proteins"

  ## remove "_" at the beginning and end
  evidence$Modified.sequence <- gsub("_", "", evidence$Modified.sequence)

  ################################################
  ### 1.2 remove the peptides including M sequence
  if(removeMpeptides){
    remove_m_sequence <- unique(evidence[grep("M", evidence$Modified.sequence), "Modified.sequence"])
    if(length(remove_m_sequence) > 0){
      evidence <- evidence[-which(evidence$Modified.sequence %in% remove_m_sequence), ]
    }
  }


  ################################################
  ## 2. remove peptides which are used in more than one protein
  ## we assume to use unique peptide
  ################################################
  if(useUniquePeptide){
    pepcount <- unique(evidence[, c("Proteins","Modified.sequence")]) ## Protein.group.IDs or Sequence
    pepcount$Modified.sequence <- factor(pepcount$Modified.sequence)

    ## count how many proteins are assigned for each peptide
    structure <- aggregate(Proteins~., data=pepcount, length)
    remove_peptide <- structure[structure$Proteins!=1, ]

    ## remove the peptides which are used in more than one protein
    if(length(remove_peptide$Proteins != 1) != 0){
      evidence <- evidence[-which(evidence$Modified.sequence %in% remove_peptide$Modified.sequence), ]
    }
  }

  ################################################
  ## 3. duplicated rows for certain feature and certain runs
  ## 	3.1) take highest intensity
  ##  3.2) take sum of intensities
  ################################################

  ## Let's find duplicates
  ## first remove NA intensity
  evidence <- evidence[!is.na(evidence$Intensity), ]


  #########################
  ### 2.1) general Label-free : one measurement for a feature and a run
  ## count the number of intensities for feature by runs

  ## take the highest intensity among duplicated or sum of intensities
  ## summaryforMultipleRows="max" or "sum
  infile_w <- .cast_maxquant_to_wide_glf(evidence, aggregateFun=summaryforMultipleRows)

  ## *** remove features which has less than 2 measurements across runs
  ## !!! for MSstats v3, we don't need to remove them.
  ## good to remove before reformatting to long-format

  if(fewMeasurements == "remove"){
    infile_w <- .remove_feature_with_few(infile_w)
  }

  ## then, go back to long-format
  # good to fill rows with NAs, then now can have balanced data-structure.
  infile_l <- .melt_maxquant_to_long_glf(infile_w)

  ## need to set 'IsotopeLabelType' because SILAC already has it.
  infile_l$IsotopeLabelType  <-  "L"


  #########################
  ### 2.2) label-free : however, several runs for a sample.

  ################################################
  ### merge all information
  colnames(infile_l)[1] <- "ProteinName"
  colnames(infile_l)[2] <- "PeptideSequence"
  colnames(infile_l)[3] <- "PrecursorCharge"

  ## Add in columns for FramentIon & ProductCharge (all values are NA)
  ## Add column for IsotopeLabelType (all "L")
  infile_l$FragmentIon <- NA
  infile_l$ProductCharge <- NA

  ## Create Condition & Bioreplicate columns; TODO: fill in with correct values
  infile_l <- merge(infile_l, annot, by=c("Raw.file", "IsotopeLabelType"))

  infile_l <- infile_l[, c(c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition", "BioReplicate", "Raw.file", "Intensity"))]
  colnames(infile_l)[9] <- "Run"

  infile_l$PeptideSequence <- factor(infile_l$PeptideSequence)
  infile_l$ProteinName <- factor(infile_l$ProteinName)

  return(infile_l)
}
