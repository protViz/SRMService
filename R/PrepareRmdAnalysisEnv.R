# this are some function to copy RMD files and preprare envs.
.scriptCopyHelperVec <- function(runscripts ){
  workdir <- getwd()
  for(scripts in runscripts){
    src_script <- file.path( find.package("SRMService") , scripts )

    dest_script <- file.path(workdir ,basename(scripts))
    message("copy ", src_script, " to ", dest_script)

    if(!file.copy(src_script , dest_script)){
      warning(paste("could not copy script file. Does file already exist?", dest_script, sep=" "))
    }
  }

  message(paste("your working directory now should contain ", length(runscripts) , "new files :\n",sep=" "))
  for(script in runscripts){
    message( basename(script) )
  }
}


.scriptCopyHelper <- function(runscript, rmfile ){
  workdir <- getwd()
  src_rmdfile <- file.path( find.package("SRMService") , rmfile )
  src_runscript <- file.path(find.package("SRMService"), runscript)
  message(src_rmdfile)
  message(src_runscript)

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
#'
#' Please see the Run_QuantTwoGroupAnalysis.R for more details
#'
#' @export
#' @examples
#' SRMService::RMD_Quant_2GrpAnalysis()
#'
RMD_Quant_2GrpAnalysis <- function(){
  .scriptCopyHelper("/samples/Run_QuantTwoGroupAnalysis.R","/reports/Grp2Analysis.Rmd" )
}

#' copies the RMD and Run files for the QuantQCReport in your working directory.
#'
#' Please see the Run_QuantQCReport.R for more details
#'
#' @export
#' @examples
#' SRMService::RMD_Quant_QCReport()
#'
RMD_Quant_QCReport <- function(){
  .scriptCopyHelper("/samples/Run_QuantQCReport.R","/reports/QCReport.Rmd" )
}

#' Copies the RMD and Run files for variable selection into your working directory.
#'
#' Please see the Run_QuantTwoGroupAnalysis.R for more details
#'
#' @export
#'
RMD_VarSelection <- function(){
  .scriptCopyHelper("/samples/Run_QuantQCReport.R","/reports/QCReport.Rmd" )
}

#' Copies the RMD and Run R file for Library generation your working directory.
#'
#' Please see the Run_specLProzor.R file in the working directory
#' for more details
#'
#' @export
#'
RMD_LibraryGen_specLProzor <- function(){
  .scriptCopyHelper("/samples/Run_specLWithProzor.R","/reports/specLWithProzor.Rmd" )
}

#' Copies the Rnw file and Run R file for old 1To1 QC into your working directory.
#'
#' Please see the Run_1To1_oldStyle.R file in your working directory,
#' for more details
#'
#' @export
#'
RMD_QC1To1_Old <- function(){
  .scriptCopyHelperVec(c("/samples/Run_1To1_oldStyle.R",
                         "/reports/MQ_sampleQC_overview.Rnw",
                         "/samples/helpers/QprotMatrixFunctions_rn_V3.R",
                         "/samples/images/LFQ_QC_workflow.pdf"))
}


