#  ZZZ.R
SIMvignette <- function()
{
	shell.exec(system.file("doc", "SIM.pdf", package="SIM"))
}

.onLoad <- function(libname, pkgname)
{	
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) 
		winMenuAddItem("Vignettes", "SIM", "SIMvignette()")
}

