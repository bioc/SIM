`heatmap.with.clac.plot` <-
function(input.region, dep.data, Normal.data, windowsize=windowsize, run.name)
{
     load(sprintf("%s/data/dep.data.only",run.name))

if(!is.vector(Normal.data)){
	msg = "Error in Normal.data, please insert a vector with either the names of the columns  or the column numbers of the normal samples in the dependent data. 
      Else, use Normal.data = FALSE"
	stop(msg)
}

DiseaseData <- dep.data.only
if(length(Normal.data) == 1) {
NormalData <- make.Normal.data(DiseaseData)
}
else{
if(is.numeric(Normal.data)) NormalData = DiseaseData[,Normal.data]
if(!is.numeric(Normal.data)) {
if(sum(match(Normal.data,colnames(dep.data),nomatch =0) == 0) > 0){
	  msg = paste("The columns of the Normal.data do not match with the column names
  	  of the dependent data, 
 	  please check the input of Normal.data")
    	  stop(msg)
 	}
      NormalData = DiseaseData[,Normal.data]

}
}

chromosome <- dget(sprintf("%s/data/dep.chr.data",run.name))
nucpos <- dget(sprintf("%s/data/dep.pos.data",run.name))

if(!file.exists(sprintf("%s/intermediate.data/regionmeans_clac",run.name)))
   {
NormalResult<-clac.preparenormal.R(DiseaseData, NormalData, Normal.Type=rep(0,ncol(NormalData)), chromosome.number=chromosome, nucleotide.position=nucpos, windowsize=windowsize, targetFDR=0.01, chromosomeOption=FALSE)
clac.result<-clac.tumorarray.R(NormalResult, tumorarrayIndex=1:ncol(DiseaseData))
dput(clac.result$RegionMean, file=sprintf("%s/intermediate.data/regionmeans_clac",run.name))
    }
RegionMean <- dget(sprintf("%s/intermediate.data/regionmeans_clac",run.name))
options(warn = -1)
sampleM <- as.matrix((RegionMean)[, 1:ncol(RegionMean)])
centro <- c(1.23e+08, 93500000, 91500000, 5.1e+07, 47700000,
            60500000, 58900000, 4.5e+07, 4.9e+07, 4e+07, 5.3e+07,
            35400000, 1.5e+07, 15600000, 1.7e+07, 3.9e+07, 2.4e+07,
            1.6e+07, 28500000, 27800000, 12200000, 11900000,
            58500000, 1e+07)
   if (substr(input.region, 1,3) == "chr")
    {
       chrom_full = strsplit(input.region, "_")[[1]][1]
       start_end = strsplit(input.region, "_")[[1]][2]
       chrom = as.integer(strsplit(chrom_full, "chr")[[1]][2])
       j = chrom
       start = as.integer(strsplit(start_end, "-")[[1]][1])
       end = as.integer(strsplit(start_end, "-")[[1]][2])
    
       abs.start = (chrom*1000000000)+start
       abs.end   = (chrom*1000000000)+end
    }
    if (!substr(input.region, 1,3) == "chr")
    {
       options(warn=-1)
       if(is.na(as.integer(input.region))) {
          abs.start = chrom.table[input.region, "Abs.start"]
          abs.end   = chrom.table[input.region, "Abs.end"]
          j <- chrom.table[input.region,"name"]
       }

       if(!is.na(as.integer(input.region))) {
          chr.pos <- chrom.table[input.region == chrom.table$name,]
          if(nrow(chr.pos)==2){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[2, "Abs.end"]
		}
		if(nrow(chr.pos)==1){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[1, "Abs.end"]
		}
          j <-  chr.pos[1,"name"]
       }
       options(warn=1)
    }
options(warn=1)
   chr <- dget(sprintf("%s/data/dep.chr.data", 
            run.name))
   nucposi <- dget(sprintf("%s/data/dep.pos.data", 
            run.name))
   jp <- 0
   abs.start.dep <- dget(sprintf("%s/data/abs.start.dep", 
            run.name))
   indices.dep = (abs.start.dep >= abs.start) & (abs.start.dep <= 
            abs.end)

   y=sampleM[indices.dep, ]
   chr <- chr[indices.dep]
   nuc <- nucposi[indices.dep]
   y[is.na(y)]=0
   posicount<-apply(y>0, 1, sum)
   negacount<-apply(y<0, 1, sum)
   cut <- 0
   M<-ncol(sampleM)+1  
   color<-rainbow(6*M)
   color[M-1] <- "yellow"
   color[M+1] <- "yellow"
   height<-M 
  xlim=range(0,max(nuc))

   nuc2 <- seq(0,length(nuc)-1,1)
       ran<-range(nuc)
      plot(rep(1,length(nuc)),nuc,type="n",xlim=c(1,-1),
    ylim=c(10,length(nuc)-11),ylab="",xlab="% gains and losses", xaxt = "n")
    for(i in 1:length(nuc2)){
      if(posicount[i]+negacount[i]>=cut)
        {
         segments(jp, nuc2[i],jp+(posicount[i])/height,nuc2[i],col=color[M-posicount[i]])

      ######if(negacount[i]-1 >=cut)
         segments(jp, nuc2[i],jp-(negacount[i])/height, nuc2[i],col=color[M+negacount[i]])
         }
         segments(jp, min(nuc2),jp, max(nuc2),col="black")
    }
   segments(jp+.5, centro[j],jp-.5,centro[j],col=6)
   segments(jp,xlim[1],jp,xlim[2])
   segments(0.5,xlim[1],0.5,xlim[2],col="gray", lty="dotted")
   segments(-0.5,xlim[1],-0.5,xlim[2],col="gray", lty="dotted")
   axis(1, c(1,0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75, -1), c(100,75,50,25,0,25,50,75,100), cex.axis = 0.75)

   
}

