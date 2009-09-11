#add colour key to plot
colourKey <- function(zlim, colRamp)
{	
    colKey <- seq(zlim[1], zlim[2], length=length(colRamp))
    image(x=colKey, z=as.matrix(colKey), col=colRamp, yaxt="n", xaxt="n", xlim=zlim, ann=FALSE)
    ax <- axTicks(1)	
    xlabMin <- substitute(xmin <= NULL, list(xmin=ax[1]))			   
    xlabMax <-  substitute(NULL >= xmax, list(xmax=ax[length(ax)]))
    xlab <- c(do.call("expression", list(xlabMin)), ax[-c(1, length(ax))], do.call("expression", list(xlabMax)))   
    axis(1, at=ax, labels=xlab)
}

#function to generate a quantsmooth plot next to the heatmap
quantsmooth.plot <- function(x, smooth.lambda, segment=min(nrow(x), 100), xlim, ylim=c(1, nrow(x)), col, ...)
{
    if(!is.numeric(x))
        stop("Data for smoothing is of wrong type!")
    
    #do quantsmooth on each sample in the dataset	
    x <- apply(x, 2, quantsmooth, smooth.lambda=smooth.lambda, segment=segment)
    
    y <- seq(nrow(x))
    
    xlim <- rev(xlim) #beware the xlim is reversed	
    
    if(missing(col))
        col <- rainbow(ncol(x))
    
    #plot the smoothing line
    matplot(x=x, y=y, col=col, type="l", lty=1, xlim=xlim, ylim=ylim,	
            xlab="log-ratio", ylab="", main="", yaxt="n", xaxt="n", yaxs="i", pch=rep(1, nrow(x)))
    axis(1, at=seq(round(xlim[1]), round(xlim[2])), labels=seq(round(xlim[1]), round(xlim[2])))	
    abline(v=0)
}

#function to produce a heatmap of the dependent features next to the zscore-heatmap
dependent.heatmap <- function(z, adjust, zlim, colRamp, ylim, ...)
{
    y <- seq_len(ncol(z))
    x <- seq_len(nrow(z))
    
    #rescale z	
    z[z > zlim[2]] <- zlim[2]
    z[z < zlim[1]] <- zlim[1]
    
    image(x, y, as.matrix(z), zlim=zlim, col=colRamp, ann=FALSE, axes=FALSE, yaxs="i", ylim=ylim)
    
    #draw the name of the samples as an label on the heatmap
    #give different colors when a subtype is given
    
    if(class(adjust) != "formula")
        stop("'adjust' is not a formula!")
    
    lwd.ticks <- par("pin")[1] * 90 / diff(par("usr")[1:2]) #90 instead of 96 in order to get a small with separation
    cex.axis <- 1
    if(length(rownames(z))> 0 )
        cex.axis <- par("mai")[1]/(max(strwidth(rownames(z), units="inches")) + 0.4)
    
    #only simple formula's
    if(length(adjust) == 2 & adjust[[2]] != 1) {
        adjust <- get(as.character(adjust[[2]]), env=attr(adjust, ".Environment"))
        adjust <- factor(adjust)
        if(nlevels(adjust) > 1) {
            #colour different levels
            for(i in 1:nlevels(adjust))
                axis(1, at=x[adjust==levels(adjust)[i]], labels=rownames(z)[adjust==levels(adjust)[i]], 
                        las=2, cex.axis=cex.axis, col=NA, col.ticks=i, lwd.ticks=lwd.ticks, lend=1)			
        }
    }
    else {
        axis(1, at=x, labels=rownames(z), las=2, col=NA, cex.axis=cex.axis, lwd.ticks=lwd.ticks, lend=1, col.ticks="deeppink4") 
    }
    
}


independent.heatmap <- function(z, zlim, sign.p.values, sign.z.scores.neg, sign.z.scores.pos, colRamp, ylab.col, 
        xlab.neg.col,	xlab.pos.col, ylim=c(1, ncol(z)), xlim=c(1, nrow(z)))
{	
    x <- do.call(":", as.list(xlim))
    y <- do.call(":", as.list(ylim))
    
    z <- z[x,]
    z <- z[,y]
    
    xlim <- c(xlim[1] - 0.5, xlim[2] + 0.5)
    ylim <- c(ylim[1] - 0.5, ylim[2] + 0.5)
    
    image(x, y, z, axes=FALSE, zlim=zlim, col=colRamp, ann=FALSE, xlim=xlim, ylim=ylim, yaxs="i")
   
    sign.p.values <- sign.p.values[y] 
    sign.z.scores.neg <- sign.z.scores.neg[x] 
    sign.z.scores.pos <- sign.z.scores.pos[x]
    
    ylab <- colnames(z)[sign.p.values]
    xlab.neg <- rownames(z)[sign.z.scores.neg]
    xlab.pos <- rownames(z)[sign.z.scores.pos]
    
    lwd.ticks <- par("pin")[2] * 90 / diff(par("usr")[3:4]) #90 instead of 96 in order to get a small with separation
          
    if(length(ylab) > 0 ) {
        ycex.axis <- par("mai")[1]/(max(strwidth(ylab, units="inches")) + 0.4)		
        axis(2, at=y[sign.p.values], labels=FALSE, col=NA, col.ticks=ylab.col, lwd.ticks=lwd.ticks, lend=1)
        axis(4, at=y[sign.p.values], labels=ylab, col=NA, col.ticks=ylab.col, lwd.ticks=lwd.ticks, las=1, cex.axis=ycex.axis, lend=1)
    }
    
    lwd.ticks <- par("pin")[1] * 90 / diff(par("usr")[1:2])
    
    if(length(xlab.neg)> 0 ) {		
        xcex.axis <- par("mai")[1]/(max(strwidth(xlab.neg , units="inches")) + 0.4)    
        axis(1, at=x[sign.z.scores.neg], labels=xlab.neg, col=NA, col.ticks=xlab.neg.col, las=2, cex.axis=xcex.axis, lwd.ticks=lwd.ticks, lend=1)
    }
    
    if(length(xlab.pos)> 0 ) {
        xcex.axis <- par("mai")[1]/(max(strwidth(xlab.pos, units="inches")) + 0.4)
        axis(1, at=x[sign.z.scores.pos], labels=xlab.pos, col=NA, col.ticks=xlab.pos.col, las=2, cex.axis=xcex.axis, lwd.ticks=lwd.ticks, lend=1)
    }    
}

heatmap <- function(dep.data, z.scores, adjusted.p.values, significance, z.threshold, colRamp, add.colRamp, add.plot, 
        adjust, add.scale, scale, smooth.lambda, input.region, ...)
{
          
    sign.p.values <- adjusted.p.values <= significance    
    
    sign.z.scores <- z.scores[sign.p.values, , drop=FALSE] 
    
    mean.z.scores <- apply(sign.z.scores, 2, mean, na.rm=TRUE)
    
    mean.z.scores[is.nan(mean.z.scores)] <- 0 #make them unsignificant
    
    sign.z.scores.neg <- mean.z.scores <= -z.threshold    
    sign.z.scores.pos <- mean.z.scores >= z.threshold
      
    #scale the zscores 
    if(missing(scale)){
        zmax <- max(0, abs(z.scores), na.rm=TRUE)
        zlim <- c(-zmax, zmax)
    }
    else		
        zlim <- sort(scale)
    
    #rescale Z-scores
    z.scores[z.scores < zlim[1]] <- zlim[1]
    z.scores[z.scores > zlim[2]] <- zlim[2]
    
    #scale the zscores 
    if(missing(add.scale)){
        xmax <- max(0, abs(as.matrix(dep.data)))
        add.scale <- c(-xmax, xmax)
    }
    else		
        add.scale <- sort(add.scale)
    
    #define colors		
    col.sign.p.values <- "gold"		
    col.sign.mean.z.scores.neg <- colRamp[1]
    col.sign.mean.z.scores.pos <- colRamp[length(colRamp)]
    
    #define layout for the different plotting options
    lts=list(smooth=list(matrix=rbind(c(1, 3), c(0, 2)), heights=c(4,1), widths=c(1,3)), 
            heatmap= list(matrix=rbind(c(1, 4), c(2, 3)), heights=c(4,1), widths=c(2,3)),
            none=list(matrix=cbind(c(2, 1)), heights=c(4,1)))
    
    #define figure margins for the different plotting options
    mars=list(smooth=list(addplot=c(5.1, 2.1, 1.1, 0.5), 
                    colKey=c(4.1, 2.5, 2.1, 6.1),  
                    heatmap=c(5.1, 0.5, 1.1, 4.1)
            ), 
            heatmap=list(addplot=c(5.1, 2.1, 1.1, 0.5), 
                    colKeyAdd=c(4.1, 4.1, 2.1, 2.5), 
                    colKey=c(4.1, 2.5, 2.1, 6.1), 
                    heatmap=c(5.1, 0.5, 1.1, 4.1)
            ), 
            none=list(colKey=c(3.1, 3.1, 2.1, 6.1), 
                    heatmap=c(5.1, 1.1, 1.1, 4.1)
            ))
    
    
    layout(lts[[add.plot]][["matrix"]], heights=lts[[add.plot]][["heights"]], widths=lts[[add.plot]][["widths"]])
    #layout.show()
    par(oma=c(0,0,3,0))
    
    ylim <- c(1 - 0.5, nrow(z.scores) + 0.5)
    arglist <- list(...)	
    if(length(arglist) != 0)
        if(any(names(arglist) %in% "ylim"))
            ylim <- arglist$ylim
    
    #add plot		
    if(add.plot == "smooth") {
        par(mar=mars[["smooth"]][["addplot"]])
        plotAdded <- quantsmooth.plot(as.matrix(dep.data), smooth.lambda=smooth.lambda, xlim=add.scale, ylim=ylim, col=add.colRamp)
    }else if(add.plot == "heatmap") {		
        par(mar=mars[["smooth"]][["addplot"]])
        plotAdded <- dependent.heatmap(t(dep.data), adjust=adjust, zlim=add.scale, colRamp=add.colRamp, ylim=ylim)
        par(mar=mars[[add.plot]][["colKeyAdd"]])
        colourKey(add.scale, add.colRamp)
    }
    
    #colour key
    par(mar=mars[[add.plot]][["colKey"]])
    colourKey(zlim, colRamp)
    
    #heatmap
    par(mar=mars[[add.plot]][["heatmap"]])
    
    independent.heatmap(z=t(z.scores), zlim=zlim, sign.p.values=sign.p.values, sign.z.scores.neg=sign.z.scores.neg,  
            sign.z.scores.pos=sign.z.scores.pos, colRamp=colRamp, ylab.col=col.sign.p.values,     
            xlab.neg.col=col.sign.mean.z.scores.neg, xlab.pos.col=col.sign.mean.z.scores.pos, ...)
    
    mtext(paste("\nHeatmap with P-value cut-off(",significance, "), Z-score threshold(", 
                    paste(-z.threshold, z.threshold, sep=", "),")\n",
                    "input region:", input.region), side=3, outer=TRUE)
    
    #still par-setting not right
    call <- match.call()
    attr(call, "env") <- parent.frame()
    invisible(call)
}

sim.plot.zoom.in <- function(call)
{
    usr <- par("usr")	
    
    cat("Click on two points in the heatmap to zoom in:\n")
    flush.console() #for windows
    
    location <- FALSE
    while(!is.list(location)) {
        location <- locator(2)
        xlim <- sort(location$x)
        ylim <- sort(location$y)
        xlim <- round(xlim) 
        ylim <- round(ylim)
        if(xlim[1] < usr[1] | xlim[2] > usr[2] | ylim[1] < usr[3] | ylim[2] > usr[4]) {			
            cat("Click inside the image!\n") 		
            location <- FALSE
        } else if(diff(xlim) == 0 | diff(ylim) == 0) {
            cat("Too small region selected!\n") 		
            location <- FALSE
        }
    }
    
    call$xlim <- xlim
    call$ylim <- ylim
    
    xlim <- c(xlim[1] - 0.5, xlim[2] + 0.5)	
    ylim <- c(ylim[1] - 0.5, ylim[2] + 0.5)	
    
    rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, border="white")	
    rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, lty=2)
    
    dev.new()
    
    invisible(eval(call, env=attr(call, "env")))	#local function environment
    
}


sim.plot.zscore.heatmap <- function(input.regions="all chrs", 
        input.region.indep=NULL, 
        method=c("full", "smooth", "window", "overlap"),
        adjust=~1, 
        significance=0.2, 
        z.threshold=3, 
        colRamp=colorRampPalette(c("red", "black", "green"))(7),
        add.colRamp=colorRampPalette(c("blue", "black", "yellow"))(7), 
        show.names.indep=FALSE, 
        show.names.dep=FALSE,
        adjust.method="BY", 
        scale , 
        add.scale,
        add.plot=c("smooth", "none", "heatmap"), 
        smooth.lambda=2,		
        pdf=TRUE,		
        run.name="analysis_results", ...)
{
    ##CHECK INPUT ARGUMENTS	              
    add.plot <- match.arg(add.plot)	
    method <- match.arg(method)
    
    #check for correct adjust	
    #check user definend z.scale c(min, max)
    
    #check significance
    if(!is.numeric(significance) | (significance < 0) | (significance > 1))
        stop("'significance' parameter should be a numeric value between 0 and 1! The input is: ", significance)
    
    ##FILES THAT ALREADY EXISTS OR WILL BE CREATED	 		
    #load
    ifDepData <- file.path(run.name, "data", "dep.data")
    ifIndepData <- file.path(run.name, "data", "indep.data")
    
    #dget
    ifIndepAbsPos <- file.path(run.name, "data", "abs.start.indep")
    ifDepAbsPos <- file.path(run.name, "data", "abs.start.dep")
    
    ifPvaluesInputRegion <- file.path(run.name, method, "gpvals.pat.c.%s")
    ifZscoresInputRegion <- file.path(run.name, method, "zmat.%s") 
    ifDepDataInputRegion <- file.path(run.name, method, "data.dep.pat.c.%s") 
    ifDepId <- file.path(run.name, "data", "dep.id.data")
    ifIndepId <- file.path(run.name, "data", "indep.id.data")
    ifDepSym <- file.path(run.name, "data", "dep.sym.data")
    ifIndepSym <- file.path(run.name, "data", "indep.sym.data")
    
    ofPDF <- file.path(run.name, "heatmap.zscores", 
            paste("Heatmap%s", method, add.plot, adjust.method, format(significance, scientific = TRUE, digits=1), ".pdf", sep="")) 
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(run.name, ifDepData, ifIndepData, ifIndepAbsPos, ifDepAbsPos,
                    ifPvaluesInputRegion, ifZscoresInputRegion, ifDepDataInputRegion))
    
    if(show.names.dep)
        checkFiles(ifDepId)
    
    if(show.names.indep)
        checkFiles(ifIndepId)
    
    cat("Producing heatmap ...\n")
    
    ##LOAD THE DATA
    #load(ifDepData)
    #load(ifIndepData)
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    abs.start.indep.whole <- dget(ifIndepAbsPos)
    abs.start.dep.whole <- dget(ifDepAbsPos)
    
    if(!is.null(input.region.indep)){
        if(length(input.region.indep) != 1)
            stop("Only one input region for the independent data is allowed. The input is now: ", paste(input.region.indep, collapse=", "))
        input.region.indep <- unlist(sapply(input.region.indep, predefinedRegions, USE.NAMES=FALSE))
        region.indep <- getGenomicRegion(input.region.indep)		
        indices.indep <- (abs.start.indep.whole >= region.indep$absolute.start) & (abs.start.indep.whole <= region.indep$absolute.end)
    }
    
    for(input.region in input.regions)
    {
        region.dep <- getGenomicRegion(input.region)
        
        cat("... input region:", input.region, "\n")
        
        #find region independent		
        indices.dep <- (abs.start.dep.whole >= region.dep$absolute.start) & (abs.start.dep.whole <= region.dep$absolute.end)
        
        #find region independent if given
        if(is.null(input.region.indep))
            indices.indep <- (abs.start.indep.whole >= region.dep$absolute.start) & (abs.start.indep.whole <= region.dep$absolute.end)
        
        if(sum(indices.dep) == 0 & sum(indices.indep) == 0) 
            stop("Either no dependent or independent features are in this region?")
        
        load(sprintf(ifZscoresInputRegion, input.region))
        z.scores <- as.matrix(z.scores)
        
        #check dimension
        if(ncol(z.scores) != sum(indices.indep))
            stop("You probably have entered a wrong 'input.region.indep' argument,", 
                    "which should be the same as in 'integrated.analysis()'.", 
                    "If you didn't enter one there, leave it out here too.\n The input is: ", input.region.indep)
        
        if(nrow(z.scores) != sum(indices.dep))
            stop("You probably have entered a wrong 'input.regions' argument,", 
                    "which should be the same as in 'integrated.analysis()'.", 
                    "\nThe input is: ", input.region)	 
        
        #get the rownames of the heatmap
        rownames(z.scores) <- 1:nrow(z.scores)
        if(show.names.dep) {						
            dep.id <- dget(ifDepId)[indices.dep]
            rownames(z.scores) <- paste(1:length(dep.id), ": ", dep.id, sep="")
            if(file.exists(ifDepSym)){				
                symb <- dget(ifDepSym)[indices.dep]
                rownames(z.scores) <- paste(1:length(dep.id), ": ", dep.id, " (", symb, ")", sep="")
            }			
        }
        
        #get the colnames of the heatmap		
        colnames(z.scores) <- 1:ncol(z.scores)
        if(show.names.indep) {						
            indep.id <- dget(ifIndepId)[indices.indep]
            colnames(z.scores) <- paste(1:length(indep.id), ": ", indep.id, sep="")
            if(file.exists(ifIndepSym)){
                symb <- dget(ifIndepSym)[indices.indep]
                colnames(z.scores) <- paste(1:length(indep.id), ": ", indep.id, " (", symb, ")", sep="")
            }			
        }
        
        p.values <- dget(sprintf(ifPvaluesInputRegion, input.region))
        
        adjusted.p.values <- p.adjust(p.values, method=adjust.method)
        
        dep.data.region <- dget(sprintf(ifDepDataInputRegion, input.region))
        
        if(pdf){						
            options(error=dev.off)            
            on.exit(options(error=NULL))
            ofPDFfull <- sprintf(ofPDF, input.region)            
            pdf(ofPDFfull, ...)
        } 
        
        h <- heatmap(dep.data.region, z.scores, adjusted.p.values, significance, z.threshold, colRamp, add.colRamp, add.plot, adjust, 
                add.scale, scale, smooth.lambda, input.region)
        
        if(pdf){                        
            cat("... heatmap plot stored with file name: \n", ofPDFfull, "\n", sep="")
            dev.off()
        }
        
    }
    invisible(h)
}


