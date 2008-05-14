`A4` <-
function(filename, landscape=T) {
  if (landscape)
    pdf(filename, width=11.69, height=8.27)
  else
    pdf(filename, width=8.27, height=11.69)
}

