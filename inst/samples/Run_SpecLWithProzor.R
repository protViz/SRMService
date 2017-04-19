library(prozor)
library(specL)
library(bibliospec)

#BLIB_REDUNDANT <- "test0.9.blib"
#BLIB_FILTERED <- "test0.9.filtered.blib"

NON_REDUNDANT <- "test0.9.filtered.blib"
REDUNDANT <- "test0.9.blib"


OUTPUTDIR = "C:/Users/wolski/prog/generateSpecLibrary/res.49562"
FASTA_FILE <- "D:/projects/p2069/data/fasta/p2069_db1_d_20160322.fasta"
MAX_IONS <- 6
MIN_IONS <- 5
MZ_ERROR <- 0.05 # e.g for Q-Exactive
MASCOT_MIN_SCORE <- 17
IRT_PEPTIDES <- rename(bibliospec::CiRTpeptides, c("iRT"="rt"))
