# this are some function to copy RMD files and preprare envs.

.scriptCopyHelper <- function(runscript, rmfile ){
  workdir <- getwd()
  src_rmdfile <- file.path( path.package("SRMService") , rmfile )
  src_runscript <- file.path(path.package("SRMService"), runscript)

  dest_rmdfile <- file.path(workdir ,basename(src_rmdfile))
  dest_runscript <- file.path(workdir, basename(src_runscript))

  if(!file.copy(src_rmdfile , dest_rmdfile)){
    warning("could not copy file. Does file already exist?")
  }else if(!file.copy(src_runscript,dest_runscript )){
    warning("could not copy file. Does file already exist?")
  }

  message("your working directory now should contain 2 files :",
          basename(src_rmdfile),
          " and ",
          basename(src_runscript),"\n" )
}

#' copies the RMD and Run files for 2 grp analysis in your working directory.
#' Please see the Run_QuantTwoGroupAnalysis.R for more details
#' @export
#' @examples
#' SRMService::RMD_Quant_2GrpAnalysis()
#'
RMD_Quant_2GrpAnalysis <- function(){
  .scriptCopyHelper("/samples/Run_QuantTwoGroupAnalysis.R","/reports/Grp2Analysis.Rmd" )
}

#' copies the RMD and Run files for the QuantQCReport in your working directory.
#' Please see the Run_QuantQCReport.R for more details
#' @export
#' @examples
#' SRMService::RMD_Quant_QCReport()
#'
RMD_Quant_QCReport <- function(){
  .scriptCopyHelper("/samples/Run_QuantQCReport.R","/reports/QCReport.Rmd" )
}

#' Copies the RMD and Run files for 2 variable selection into your working directory.
#' Please see the Run_QuantTwoGroupAnalysis.R for more details
#' @export
#'
RMD_VarSelection <- function(){
  .scriptCopyHelper("/samples/Run_QuantQCReport.R","/reports/QCReport.Rmd" )
}

#' Copies the RMD and Run files for 2 grp analysis in your working directory.
#' Please see the Run_QuantTwoGroupAnalysis.R for more details
#'
RMD_LibraryGen_specLProzor <- function(){
  .scriptCopyHelper("/samples/Run_specLProzor.R","/reports/specLWithProzor.Rmd" )
}

