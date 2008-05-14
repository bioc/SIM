`create.color.vector` <-
function(my.class, levels=NULL, color.vector=NULL, my.start=0, my.end=0.7) {
  
  # Coerce "my.class" into a factor.
  my.class = factor(my.class, levels)
  
  # Create a list that contains as many colors as "my.class" has levels.
  # "start" and "end" shift the hue if necessary.
  if ( is.null(color.vector) ) 
    color.vector = rainbow(nlevels(my.class), start=my.start, end=my.end)  
  
  # Create a list with as many entries as "my.class".
  my.class.color = rep(0, length(my.class))
  
  # Set each entry of "my.class.color" to the correct color in "color.vector".
  for(i in 1:length(my.class))
    my.class.color[i] = color.vector[ my.class[i] ]
  
  
  # Return the vector with colors whose entries correspond to the passed
  # argument.
  return( list(my.class.color, color.vector) )
}

