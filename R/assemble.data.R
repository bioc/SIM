`assemble.data` <-
function (dep.data = acgh.data, indep.data = expr.data, ann.dep = colnames(acgh.data)[1:4], 
    ann.indep = colnames(expr.data)[1:4], dep.id = "ID", dep.chr = "CHROMOSOME", 
    dep.pos = "STARTPOS", dep.symb=FALSE, indep.id = "ID", indep.chr = "CHROMOSOME", 
    indep.pos = "STARTPOS", indep.symb=FALSE, overwrite = FALSE, run.name = NULL) {
	Abs.start = FALSE
    if (is.null(run.name)) {
        run.name = "analysis_results"
    }
    if (!is.data.frame(dep.data) | !is.data.frame(indep.data)) {
        msg = paste("The dependent and independent datasets should be data.frames.")
        stop(msg)
    }
    if (file.exists(run.name) & !overwrite) {
        msg = paste("The directory", run.name, "should not exist or the argument 'overwrite' should be \n   \t  set to TRUE.")
        stop(msg)
    }
    if (!length(dep.id) == 1 | !length(indep.id) == 1 | !length(dep.pos) ==
        1 | !length(dep.symb) == 1 | !length(indep.pos) == 1 |
        !length(dep.chr) == 1 | !length(indep.chr) == 1 | !length(indep.symb) ==
        1 | !length(Abs.start) == 1) {
        msg = paste("the inserted names of the ID's, positions, chromosomes and Abs.start \n\t\t should be single strings, not a list of strings or being empty")
        stop(msg)
    }
    sprintf("run.name %s", run.name)
    options(warn = -1)
    dir.create(run.name)
    dir.create(sprintf("%s/data/", run.name))
    dir.create(sprintf("%s/pvalue_plots/", run.name))
    dir.create(sprintf("%s/intermediate.data/", run.name))
    dir.create(sprintf("%s/heatmap_zscores/", run.name))
    dir.create(sprintf("%s/top.indep.features/", run.name))
    dir.create(sprintf("%s/top.dep.features/", run.name))
    options(warn = 1)
    if (!dep.pos %in% colnames(dep.data) | !indep.pos %in% colnames(indep.data)) {
        msg = paste("The dep.pos and/or indep.pos does not match with the column names\n \t     of the dependent and/or independent data, \n  \t    please check the names of the dep.pos and indep.pos")
        stop(msg)
    }
    dep.pos.data = dep.data[, dep.pos]
    indep.pos.data = indep.data[, indep.pos]
    if (!dep.chr %in% colnames(dep.data) | !indep.chr %in% colnames(indep.data)) {
        msg = paste("The dep.chr and/or indep.chr does not match with the column names\n    \t  of the dependent and/or independent data, \n   \t  please check the names of the dep.chr and indep.chr")
        stop(msg)
    }
    dep.chr.data = dep.data[, dep.chr]
    indep.chr.data = indep.data[, indep.chr]
    if (dep.symb != FALSE) {
        if (!dep.symb %in% colnames(dep.data)) {
            msg = paste("The dep.symb does not match with the column names\n    \t  of the dependent , \n   \t  please check the names of the dep.symb")
            stop(msg)
        }
        dep.symb.data = dep.data[, dep.symb]
        
    }
      if (indep.symb != FALSE) {
          if (!indep.symb %in% colnames(indep.data))
          {
             msg = paste("The indep.symb does not match with the column names\n    \t  of the independent , \n   \t  please check the names of the indep.symb")
            stop(msg)
            }
            indep.symb.data = indep.data[, indep.symb]
            }
    
    if (!Abs.start == FALSE) {
        if (!Abs.start %in% colnames(dep.data) | !Abs.start %in%
            colnames(indep.data)) {
            msg = paste("The Abs.start does not match with the column names\n      \t\t\t of the dependent and/or independent data, \n      \t\t\t please check the name of the Abs.start")
            stop(msg)
        }
        Abs.start.indep = indep.data[, Abs.start]
        Abs.start.dep = dep.data[, Abs.start]
    }
    if (Abs.start == FALSE) {
        Abs.start.indep = get.abs.start(indep.chr.data, indep.pos.data,
            "independent")
        writeln()
        Abs.start.dep = get.abs.start(dep.chr.data, dep.pos.data,
            "dependent")
    }
    Abs.start.indep2 <- Abs.start.indep[order(Abs.start.indep)]
    Abs.start.dep2 <- Abs.start.dep[order(Abs.start.dep)]
    dput(Abs.start.indep2, file = sprintf("%s/data/abs.start.indep",
        run.name))
    dput(Abs.start.dep2, file = sprintf("%s/data/abs.start.dep",
        run.name))
    dep.pos.data <- dep.pos.data[order(Abs.start.dep)]
    indep.pos.data <- indep.pos.data[order(Abs.start.indep)]
    dput(dep.pos.data, file = sprintf("%s/data/dep.pos.data",
        run.name))
    dput(indep.pos.data, file = sprintf("%s/data/indep.pos.data",
        run.name))
    dep.chr.data <- dep.chr.data[order(Abs.start.dep)]
    indep.chr.data <- indep.chr.data[order(Abs.start.indep)]
    dput(dep.chr.data, file = sprintf("%s/data/dep.chr.data",
        run.name))
    dput(indep.chr.data, file = sprintf("%s/data/indep.chr.data",
        run.name))
    if (dep.symb != FALSE) {
        dep.symb.data = dep.symb.data[order(Abs.start.dep)]
        dput(dep.symb.data, file = sprintf("%s/data/dep.symb.data",
        run.name))
    }
    if (indep.symb != FALSE) {
       indep.symb.data = indep.symb.data[order(Abs.start.indep)]
       dput(indep.symb.data, file = sprintf("%s/data/indep.symb.data",
       run.name))
       }
    dep.data <- dep.data[order(Abs.start.dep), ]
    indep.data <- indep.data[order(Abs.start.indep), ]
    save(dep.data, file = sprintf("%s/data/dep.data", run.name))
    save(indep.data, file = sprintf("%s/data/indep.data", run.name))
    if (is.integer(ann.indep))
        ann.indep <- colnames(indep.data)[ann.indep]
    if (is.integer(ann.dep))
        ann.dep <- colnames(dep.data)[ann.dep]
    if (sum(match(ann.dep, colnames(dep.data), nomatch = 0) ==
        0) > 0 | sum(match(ann.indep, colnames(indep.data), nomatch = 0) ==
        0) > 0) {
        msg = paste("The annotation names of the independent and/or dependent data do not match \t\t  with the column names\n  \t  of the dependent and/or independent data, \n \t  please check the names of the ann.dep and ann.indep")
        stop(msg)
    }
    ann.dep = dep.data[, ann.dep]
    ann.indep = indep.data[, ann.indep]
    dput(ann.dep, file = sprintf("%s/data/ann.dep.data", run.name))
    dput(ann.indep, file = sprintf("%s/data/ann.indep.data",
        run.name))
    if (!dep.id %in% colnames(dep.data) | !indep.id %in% colnames(indep.data)) {
        msg = paste("The dep.id and/or indep.id does not match with the column names\n     \t    of the dependent and/or independent data, \n      \t    please check the names of the dep.id and indep.id")
        stop(msg)
    }
    dep.id.data = dep.data[, dep.id]
    indep.id.data = indep.data[, indep.id]
    dep.id.data[dep.id.data == " "] <- NA
    indep.id.data[indep.id.data == " "] <- NA
    dput(dep.id.data, file = sprintf("%s/data/dep.id.data", run.name))
    dput(indep.id.data, file = sprintf("%s/data/indep.id.data",
        run.name))
}

