`collapse.bands` <-
function(chr.bands, level) {

  # Remove the column "stain" if present.
  idx = match("stain", colnames(chr.bands))
  
  if ( !is.na(idx) )
    chr.bands = chr.bands[, -idx]
  
  # Add the columns "Abs.start" and "Abs.end" if necessary.
  idx = match("Abs.start", colnames(chr.bands))
  
  if ( is.na(idx) ) {
    col.names = c(colnames(chr.bands), "Abs.start", "Abs.end")
    
    abs.start = chr.bands$seq_region_start + as.integer(chr.bands$name) * 10^9
    abs.end   = chr.bands$seq_region_end   + as.integer(chr.bands$name) * 10^9
  
    chr.bands = cbind(chr.bands, abs.start, abs.end)
    colnames(chr.bands) = col.names
    
  } 

  
  # Trim the band. 0 is a special case, we'll only return a single row for
  # the entire genome.
  if (level == 0) {
    chrn  = 0
    start = min(chr.bands$Abs.start)
    end   = max(chr.bands$Abs.end)
    band  = ''
    
    ret = data.frame(chrn, 1, end, band, start, end)
    colnames(ret) = colnames(chr.bands)
    
    return(ret) 
  }
  
  
  chr.bands$band = substr(chr.bands$band, 1, level-1) 
  chromosomes = factor(chr.bands$name)
  chr.bands.new = data.frame()
  
  for (chrn in levels(chromosomes)) {
    bands = factor(chr.bands[chr.bands$name == chrn, "band"])
    
    for (band in levels(bands)) {     
      indices = (chr.bands$name == chrn) & (chr.bands$band == band)
      
      start = min(chr.bands[indices, "seq_region_start"])
      end   = max(chr.bands[indices, "seq_region_end"])
            
      abs.start = min(chr.bands[indices, "Abs.start"])
      abs.end   = max(chr.bands[indices, "Abs.end"])
      
      row = data.frame(chrn, start, end, band, abs.start, abs.end)
      chr.bands.new = rbind(chr.bands.new, row)
    }
  }   
  
  colnames(chr.bands.new) = colnames(chr.bands)
  rownames(chr.bands.new) = paste(chr.bands.new$name, chr.bands.new$band, sep='')
  chr.bands.new = chr.bands.new[order(chr.bands.new$Abs.start), ]
  
  return(chr.bands.new)
}

