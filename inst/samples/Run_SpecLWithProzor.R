#
# If you are going to use results produced by the scripts please do cite the
# SRMSerivce R package by providing the following reference
# www.github.com/protViz/SRMService
# by W.E. Wolski, J. Grossmann, C. Panse
#

library(prozor)
library(specL)
library(bibliospec)

# you start with a folder containing your dat file.


# find out which fasta file was used. You will need this file to annotate the peptides
DATFILE <- "F251145.dat"
xx <-readLines(DATFILE, n=1000)
grep("fastafile=",xx,value = TRUE)

# path to fasta file
FASTA_FILE <- "D:/projects/p2244_MilenaS_PN/library/fgcz_9606_d_reviewed_cnl_contaminantNoHumanCont_20161209.fasta"


# How your library should be called and where stored:
OUTPUTDIR = "output"
NON_REDUNDANT <- "F251145.filtered.sqlite3"
REDUNDANT <- "F251145.redundant.sqlite3"
SWATH_LIBRARY <- "test_MeOHMerged.tsv"

# ion library parameters
MAX_IONS <- 6
MIN_IONS <- 5
MZ_ERROR <- 0.05 # e.g for Q-Exactive
MASCOT_MIN_SCORE <- 17
IRT_PEPTIDES <- plyr::rename(bibliospec::CiRTpeptides, c("iRT"="rt"))


# generate bibliospec files
bb<-bibliospec::Blib()
bb$build(idfiles = DATFILE, outfile = REDUNDANT,cutoff = 0.9)
bb$filter(infile = REDUNDANT, outfile = NON_REDUNDANT)


rmarkdown::render("specLWithProzor.Rmd")


