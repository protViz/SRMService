#R

.onAttach <- function(lib, pkg){
	if(interactive()){
		version <- packageVersion('SRMservice')
		packageStartupMessage("Package 'SRMservice' version ", version)
	  invisible()
	}
}
