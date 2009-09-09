#  ZZZ.R
.onAttach <- function(lib, pkg)
{	
    #load compiled C-code	
    library.dynam("SIM", pkg, lib)
    
    #for using the genes x samples format of the data	
    gt.options(transpose=TRUE)
    
    #add vignette to windowsmenu
    #if(interactive() && .Platform$OS.type == "windows" &&
    #        .Platform$GUI == "Rgui"){
    #    addVigs2WinMenu("SIM")
    #}
}

