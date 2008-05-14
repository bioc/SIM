`acgh.heatmap.with.heatmap` <-
function(subtype, input.region, scale.heatmap , run.name)
{
load(sprintf("%s/data/dep.data.only",run.name))
input.region.name <- input.region

    if (substr(input.region.name, 1,3) == "chr")
    {
       chrom_full = strsplit(input.region.name, "_")[[1]][1]
       start_end = strsplit(input.region.name, "_")[[1]][2]
       chrom = as.integer(strsplit(chrom_full, "chr")[[1]][2])
    
       start = as.integer(strsplit(start_end, "-")[[1]][1])
       end = as.integer(strsplit(start_end, "-")[[1]][2])
    
       abs.start = (chrom*1000000000)+start
       abs.end   = (chrom*1000000000)+end
    }
    else if (!substr(input.region.name, 1,3) == "chr")
    {
       options(warn=-1)
       if(is.na(as.integer(input.region.name))) {
          abs.start = chrom.table[input.region.name, "Abs.start"]
          abs.end   = chrom.table[input.region.name, "Abs.end"]
       }

       if(!is.na(as.integer(input.region.name))) {
          chr.pos <- chrom.table[input.region.name == chrom.table$name,]

          if(nrow(chr.pos)==2){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[2, "Abs.end"]
		}
		if(nrow(chr.pos)==1){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[1, "Abs.end"]
		}
       }
       options(warn=1)
    }
    
  abs.start.dep <- dget(sprintf("%s/data/abs.start.dep", 
            run.name))
        indices.dep = (abs.start.dep >= abs.start) & (abs.start.dep <= abs.end)

       if (sum(indices.dep) > 2) {
    
      acgh.data.only = as.matrix(dep.data.only[indices.dep, ])
      
     }

x = 1:ncol(acgh.data.only)
y = 1:nrow(acgh.data.only)
z = t(acgh.data.only)
  # this will transpose the graph. it's used instead of the normal
    # x=1:nrow(acgh.data.only); y=1:ncol(acgh.data.only); z=acgh.data.only

if(length(scale.heatmap) == 1) 
{
zmax = max(abs(acgh.data.only))
scale.heatmap = c(-zmax, zmax)
}

# colors
colors = maPalette(low="red", mid="black", high="green")

for(rows in 1:nrow(z)){
z[rows,(z[rows,] > scale.heatmap[2])] = scale.heatmap[2]
z[rows,(z[rows,] < scale.heatmap[1])] = scale.heatmap[1] }

image(x,y,z,zlim=scale.heatmap,col=colors,xaxt='n',yaxt='n',xlab='',ylab='')
samples <- colnames(acgh.data.only)

if(length(samples) > 30) cex.ax = 0.5
if(length(samples) <= 30) cex.ax = 1

# label the columns. for example for chromosome 1:
axis(1,at=1:length(samples),labels=samples,las=2, cex.axis = cex.ax)

if(length(subtype) == 1)
{
if(subtype != FALSE){
 stop("insert a correct subtype, either a vector with the same length as 'samples' inserted in 'integrated.analysis' or FALSE")
}
}

if(length(subtype) == length(samples))
{
subtype <- factor(subtype)
for(i in 1:nlevels(subtype))
{
sub = grep(levels(subtype)[i], subtype, T)
axis(1, at=sub,labels=samples[subtype == levels(subtype)[i]], las=2, col.axis=i+1, cex.axis = cex.ax)
}
}
}

