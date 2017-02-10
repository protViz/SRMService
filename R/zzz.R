#R

.onAttach <- function(lib, pkg){
	if(interactive()){
		version <- packageVersion('SRMService')
		packageStartupMessage("Package 'SRMService' version ", version)
	  invisible()
	}
}
