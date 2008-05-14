`integrated.analysis` <-
function (samples, input.regions = "all chrs", 
    adjust = FALSE, zscores = FALSE, method = c("auto", "asymptotic", 
        "permutations", "gamma"), run.name = NULL) 
{
    data(chrom.table)
    miss_arms = c("13p", "14p", "15p", "21p", "22p")
    writeln(sprintf("Running the integrated analysis for: %s", 
        input.regions))
    options(warn = -1)
    if (input.regions == "all chrs") {
        input.regions = c(1:24)
    }
    else if (input.regions == "all chrs auto") {
        input.regions = c(1:22)
    }
    else if (input.regions == "all arms") {
        input.regions = rownames(chrom.table)
        input.regions <- input.regions[!(input.regions %in% miss_arms)]
    }
    else if (input.regions == "all arms auto") {
        input.regions = rownames(chrom.table)[1:44]
        input.regions <- input.regions[!(input.regions %in% miss_arms)]
    }
    options(warn = 1)
    if (is.null(run.name)) {
        run.name = "analysis_results"
    }
    if (!file.exists(run.name)) {
        msg = paste("The directory", run.name, "does not exist, please insert a correct run.name.")
        stop(msg)
    }
    if (!file.exists(sprintf("%s/data/dep.data", run.name))) {
        msg = "please first run the read.data() function to read in the data necessary for this function"
        stop(msg)
    }
    load(sprintf("%s/data/dep.data", run.name))
    load(sprintf("%s/data/indep.data", run.name))
    if (is.integer(samples)) {
        samples.indep <- colnames(indep.data)[samples]
        samples.dep <- colnames(dep.data)[samples]
        if (sum(match(samples.dep, colnames(dep.data), nomatch = 0) == 
            0) > 0 | sum(match(samples.indep, colnames(indep.data), 
            nomatch = 0) == 0) > 0) {
            stop("The inserted samples do not match with the column names of the dependent and/or  \tindependent data, please check samples input")
        }
        indep.data.only = indep.data[, samples.indep]
        dep.data.only = dep.data[, samples.dep]
    }
    if (!is.integer(samples)) {
        if (sum(match(samples, colnames(dep.data), nomatch = 0) == 
            0) > 0 | sum(match(samples, colnames(indep.data), 
            nomatch = 0) == 0) > 0) {
            stop("The inserted samples do not match with the column names of the dependent and/or  \tindependent data, please check samples input2")
        }
        indep.data.only = indep.data[, samples]
        dep.data.only = dep.data[, samples]
    }
    save(dep.data.only, file = sprintf("%s/data/dep.data.only", run.name))
    save(indep.data.only, file = sprintf("%s/data/indep.data.only", run.name))
    if (sum(is.na(dep.data.only)) > 0) {
        stop("Missing values (NA) are not allowed in the dep.data, please run the function\n\timpute.nas.by.surrounding(), to replace the NA's by the median of the surrounding features")
    }
    if (sum(match(colnames(indep.data.only), colnames(dep.data.only), 
        nomatch = 0) == 0) > 0 | sum(match(colnames(dep.data.only), 
        colnames(indep.data.only), nomatch = 0) == 0) > 0) {
        msg = paste("The column names of the independent and dependent data do not match, \n      please check the column names")
        stop(msg)
    }
    for (row in input.regions) {
        writeln(sprintf("Running for: %s", row))
        if (row == "whole genome auto") {
            abs.start = min(chrom.table[1:44, ]$Abs.start)
            abs.end = max(chrom.table[1:44, ]$Abs.end)
        }
        else if (row == "whole genome") {
            abs.start = min(chrom.table$Abs.start)
            abs.end = max(chrom.table$Abs.end)
        }
	else if (substr(row, 1, 3) == "chr"){
	if(!(regexpr("_", row) == 5 | regexpr("_", row) == 6) |  regexpr("-", row) < 0) {
		stop(sprintf("the inserted input.regions %s is not correct, please insert a correct input region", row))
		}
		chrom_full = strsplit(row, "_")[[1]][1]
            start_end = strsplit(row, "_")[[1]][2]
            chrom = as.integer(strsplit(chrom_full, "chr")[[1]][2])
            start = as.integer(strsplit(start_end, "-")[[1]][1])
            end = as.integer(strsplit(start_end, "-")[[1]][2])
            abs.start = (chrom * 1e+09) + start
            abs.end = (chrom * 1e+09) + end
	}
	else if (!substr(row, 1, 3) == "chr"){
	    if(!nchar(row) < 4 | !regexpr("[0-9]", row)==1 | !(regexpr("[pq]", row) == 2 | regexpr("[pq]", row) == 3 | regexpr("[a-z]", row) < 0)) {
	     stop(sprintf("the inserted input.regions %s is not correct, please insert a correct input region", row))
		}
             options(warn = -1)
            if (is.na(as.integer(row))) {
                abs.start = chrom.table[row, "Abs.start"]
                abs.end = chrom.table[row, "Abs.end"]
            }
            if (!is.na(as.integer(row))) {
                chr.pos <- chrom.table[row == chrom.table$name, 
                  ]
                if (nrow(chr.pos) == 2) {
                  abs.start = chr.pos[1, "Abs.start"]
                  abs.end = chr.pos[2, "Abs.end"]
                }
                if (nrow(chr.pos) == 1) {
                  abs.start = chr.pos[1, "Abs.start"]
                  abs.end = chr.pos[1, "Abs.end"]
                }
            }
            options(warn = 1)
        }
        abs.start.indep <- dget(sprintf("%s/data/abs.start.indep", 
            run.name))
        abs.start.dep <- dget(sprintf("%s/data/abs.start.dep", 
            run.name))
        indices.indep = (abs.start.indep >= abs.start) & (abs.start.indep <= 
            abs.end)
        indices.dep = (abs.start.dep >= abs.start) & (abs.start.dep <= 
            abs.end)
        if (!sum(indices.indep) > 0 | !sum(indices.dep) > 0) {
            msg = paste("There are no no features to run for, please check your data or input regions")
            stop(msg)
        }
        if (sum(indices.indep) > 0 | sum(indices.dep) > 0) {
            indep.matrix = as.matrix(indep.data.only[indices.indep, 
                ])
            dep.matrix = as.matrix(dep.data.only[indices.dep, 
                ])
            dep.pos.data = dget(sprintf("%s/data/dep.pos.data", 
                run.name))
            start.input.region = dep.pos.data[indices.dep]
            my = run.ia(indep.matrix, dep.matrix, adjust, zscores, 
                method = method)
            dput(start.input.region, file = sprintf("%s/intermediate.data/start.input.region.%s", 
                run.name, row))
            dput(dep.matrix, file = sprintf("%s/intermediate.data/data.dep.pat.c.%s", 
                run.name, row))
            dput(indep.matrix, file = sprintf("%s/intermediate.data/data.indep.pat.c.%s", 
                run.name, row))
            if (!zscores) {
                dput(my[[1]], file = sprintf("%s/intermediate.data/gpvals.pat.c.%s", 
                  run.name, row))
            }
            if (zscores) {
                dput(my[[1]], file = sprintf("%s/intermediate.data/gpvals.pat.c.%s", 
                  run.name, row))
                dput(my[[2]], file = sprintf("%s/intermediate.data/my.zmat.%s", 
                  run.name, row))
                dput(my[[3]], file = sprintf("%s/intermediate.data/zmat.%s", 
                  run.name, row))
                dput(my[[4]], file = sprintf("%s/intermediate.data/zmat.col.%s", 
                  run.name, row))
            }
        }
    }
}

