`writeln` <-
function(txt="", ...) {

  if (length(list(...)) > 0)
    txt = sprintf(txt, ...)
  
  write(txt, file=stdout())
  flush.console()
}

