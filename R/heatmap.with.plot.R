`heatmap.with.plot` <-
function(dep.data, input.region = input.region, plot.method = plot.method,Normal.data = Normal.data, subtype = subtype, acgh.heatmap.scale = acgh.heatmap.scale,run.name = run.name, lambda = lambda, windowsize = windowsize, x=x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
    distfun = dist, hclustfun = hclust, reorderfun = function(d,
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv,
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,
    margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 +
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,
    verbose = getOption("verbose"), ...)
{
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv))
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv))
        Colv <- colMeans(x, na.rm = na.rm)
    rowInd <- 1:nr

    colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow))
        if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol))
        if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(1, 4)
    lhei <- c((0.05) + if (!is.null(main)) 0.2 else 0,
        4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) !=
            nc)
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) !=
            nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei,
            "; lmat=\n")
        
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none")
        x <- t(x)
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    plot.method = match.arg(plot.method, c("smooth", "clac", "none", "heatmap"))
    if (plot.method == "smooth")
        quantsmooth.heatmap(input.region, dep.data, lambda = lambda, run.name) 
    if (plot.method == "clac")
        heatmap.with.clac.plot(input.region, dep.data,Normal.data,windowsize = windowsize, run.name)
    if (plot.method == "none")
        frame() 
    if(plot.method == "heatmap")
	 acgh.heatmap.with.heatmap(input.region, subtype = subtype, acgh.heatmap.scale, run.name)
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (!is.null(main)){
	 if(nchar(main) < 5)
        title(main, cex.main = 1.5 * op[["cex.main"]], col.main = "blue")
	else title(paste("                          ",main), cex.main = 1.25, col.main= "blue")
}
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

