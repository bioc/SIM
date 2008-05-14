`get.indep.zscores` <-
function(arg,run.name) {

  filename = sprintf("%s/intermediate.data/zmat.%s",run.name, arg)
  return(dget(filename))
}

