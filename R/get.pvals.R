`get.pvals` <-
function(input.regions,run.name) {
filename = sprintf("%s/intermediate.data/gpvals.pat.c.%s", run.name,input.regions)
return( dget(filename) );
}

