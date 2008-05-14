`sim.plot.pvals.on.genome` <-
function (input.regions = "all chrs", adjust.method = c("BY", 
    "BH", "raw"), pdf = TRUE, run.name = NULL,...) 
{
while(names(dev.cur()) == "pdf") dev.off()
    data(chrom.table)
    miss_arms = c("13p", "14p", "15p", "21p", "22p")
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
    chr.bands = chrom.table
    chr.bands.split = split(chr.bands$seq_region_end, factor(chr.bands$name))
    lens = sapply(chr.bands.split, max)
    nchrs = length(lens)
    maxdwidth = max(lens)
    filename = sprintf("%s/whole_genome_plot.pdf", run.name)
    if(pdf==TRUE)
      pdf(filename, width = 11.69, height = 8.27)
    plot(c(0, maxdwidth), c(0.5, nchrs + 0.5), type = "n", ylab = "Chromosome", 
        xlab = "", axes = FALSE, las = 2, ...)
    title("Significantly associated features")
    def.legend = c("p > 0.20", "0.05 < p <= 0.20", "p <= 0.05")
    def.fill = c("#9f9f9f", "#0088ff", "#00ffff")
    legend("bottomright", legend = def.legend, cex = 0.5, fill = def.fill)
    axis(2, 1:nchrs, nchrs:1, las = 2)
    for (i in 1:nchrs) lines(c(0, lens[i]), c(nchrs + 1 - i, 
        nchrs + 1 - i), lwd = 1)
    arms = collapse.bands(chr.bands, 2)
    p.arms = arms[arms$band == "p", ]
    for (idx in 1:nrow(p.arms)) points(p.arms[idx, "seq_region_end"], 
        nchrs + 1 - idx, col = "purple", bg = "purple", pch = 19)
    hheight = 0.4
    for (input.region.name in input.regions) {
        if (substr(input.region.name, 1, 3) == "chr") {
            chrom_full = strsplit(input.region.name, "_")[[1]][1]
            start_end = strsplit(input.region.name, "_")[[1]][2]
            chrom = as.integer(strsplit(chrom_full, "chr")[[1]][2])
            chr.pos <- chrom.table[chrom == chrom.table$name, 
                ]
            start = as.integer(strsplit(start_end, "-")[[1]][1])
            end = as.integer(strsplit(start_end, "-")[[1]][2])
		x0 = start
            x1 = end 
            abs.start = (chrom * 1e+09) + start
            abs.end = (chrom * 1e+09) + end
        }
        else if (!substr(input.region.name, 1, 3) == "chr") {
            options(warn = -1)
            if (is.na(as.integer(input.region.name))) {
                abs.start = chrom.table[input.region.name, "Abs.start"]
                abs.end = chrom.table[input.region.name, "Abs.end"]
                x0 = chrom.table[input.region.name, "seq_region_start"]
                x1 = chrom.table[input.region.name, "seq_region_end"]
            }
            if (!is.na(as.integer(input.region.name))) {
                chr.pos <- chrom.table[input.region.name == chrom.table$name, 
                  ]
		if(nrow(chr.pos)==2){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[2, "Abs.end"]
                   x0 = chr.pos[1, "seq_region_start"]
                   x1 = chr.pos[2, "seq_region_end"]
		}
		if(nrow(chr.pos)==1){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[1, "Abs.end"]
                   x0 = chr.pos[1, "seq_region_start"]
                   x1 = chr.pos[1, "seq_region_end"]
		}
            }
            options(warn = 1)
        }
        y0 = nchrs + 1 - floor(abs.start/10^9)
        y1 = nchrs + 1 - floor(abs.end/10^9)
        points(x0, y0, col = "orange", pch = 24, bg = "orange", 
            cex = 0.7)
        points(x1, y1, col = "orange", pch = 25, bg = "orange", 
            cex = 0.7)
        abs.start.dep <- dget(sprintf("%s/data/abs.start.dep", 
            run.name))
        indices.dep = (abs.start.dep >= abs.start) & (abs.start.dep <= 
            abs.end)
        if (sum(indices.dep)) {
	    dep.pos.data <- dget(sprintf("%s/data/dep.pos.data", run.name))
	    dep.chr.data <- dget(sprintf("%s/data/dep.chr.data", run.name))
            x0 = x1 = dep.pos.data[indices.dep]
            y0 = nchrs + 1 - dep.chr.data[indices.dep] - 
                hheight
            y1 = nchrs + 1 - dep.chr.data[indices.dep] + 
                hheight

            segments(x0, y0, x1, y1, col = "#9f9f9f")
            raw.pvals = try(get.pvals(input.region.name, run.name = run.name), 
                silent = T)
            if (class(raw.pvals) != "try-error") {
                adjust.method = match.arg(adjust.method, c("BY", 
                  "BH", "raw"))
                if (adjust.method == "BY") 
                  p.values = adj.byfdr(raw.pvals)
                else if (adjust.method == "BH") 
                  p.values = adj.bhfdr(raw.pvals)
                else if (adjust.method == "raw") 
                  p.values = raw.pvals
                summary(p.values)
                indices = (p.values < 0.2)
                if (sum(indices) > 0) 
                  segments(x0[indices], y0[indices], x1[indices], 
                    y1[indices], col = "#0088ff")
                dget
                indices = (p.values < 0.05)
                if (sum(indices) > 0) 
                  segments(x0[indices], y0[indices], x1[indices], 
                    y1[indices], col = "#00ffff")
            }
        }
    }
    if(pdf==TRUE){
      dev.off()
      writeln()
      writeln("**************************************************")
      writeln("The results are stored in the following file:")
      writeln(filename)
      writeln("**************************************************")
      writeln()
    }
}

