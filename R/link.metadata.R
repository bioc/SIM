`link.metadata` <-
function(data = expr.data, col.ID.link = 1, chr = as.list(hgu133plus2CHR), chrloc = as.list(hgu133plus2CHRLOC), symbol = as.list(hgu133plus2SYMBOL))
{
exp.chrloc <- chrloc[!is.na(chrloc)] 
exp.chr <- chr[!is.na(chr)] 
exp.symb <- symbol[!is.na(symbol)] 

exp.chrloc <- exp.chrloc[names(exp.chrloc) %in% data[,col.ID.link]]
exp.chr <- exp.chr[names(exp.chr) %in% data[,col.ID.link]]
exp.symb <- exp.symb[names(exp.symb) %in% data[,col.ID.link]]
exp.chr <- exp.chr[names(exp.chr) %in% names(exp.chrloc)]
exp.chr <- exp.chr[names(exp.chr) %in% names(exp.symb)]

if(length(exp.chrloc) ==0 | length(exp.chr) == 0 | length(exp.symb) == 0)
{
stop("there are no similarities between the inserted chr and col.ID.link or chrloc and col.ID.link or symbol and col.ID.link")
}
exp.chr3 <- unlist(exp.chr,use.names=T)

extra.ids <- names(exp.chr3)[match(names(exp.chr3),names(exp.chr),nomatch=0)==0]
exp.chr4 <-  exp.chr3[match(names(exp.chr3),extra.ids,nomatch=0)==0]

exp.chr4 <- exp.chr4[order(names(exp.chr4))]

if(!nlevels(factor(names(exp.chr4))) == length(exp.chr4))
{
exp.chr5 <- aggregate(exp.chr4 , list(
Region=names(exp.chr4)), mean)
}

exp.chrloc3 <- unlist(exp.chrloc,use.names=T)

# so that we have no correspondence left
my.names <- rep("",length(exp.chrloc3))

for(xi in 1:length(exp.chrloc3))
 {
  my.names[xi] <- strsplit(names(exp.chrloc3)[xi],"\\.")[[1]][1]
 }
names(exp.chrloc3) <- my.names

exp.chrloc3 <- exp.chrloc3[order(names(exp.chrloc3))]

exp.chrloc4 <- aggregate(exp.chrloc3, list(
Region=names(exp.chrloc3)), mean)

exp.symb2 <- unlist(exp.symb,use.names=T) 


exp.chrloc4[,2] <- abs(exp.chrloc4[,2])
# getting only probes present in both lists
my.names <- exp.chrloc4[,1]
common.ann <- intersect(names(exp.chr4),my.names)
sel.chr <- match(names(exp.chr4),common.ann,nomatch=0)
sel.chrloc <- match(levels(factor(my.names)),common.ann,nomatch=0)
exp.chrloc <- exp.chrloc4[sel.chrloc>0,]
exp.chrloc.names <- levels(factor(my.names))[sel.chrloc>0]
exp.chr <- exp.chr4[sel.chr>0]
chr.exp <- exp.chr[order(names(exp.chr))]

expr_inc_chrloc <- merge(exp.chrloc, data,by.x=1,by.y=col.ID.link)

exp.chr.frame <- data.frame(names(exp.chr),exp.chr)
exp.chr.frame <- exp.chr.frame[order(exp.chr.frame[,1]),]

expr_inc_chr <- merge(exp.chr.frame, expr_inc_chrloc,by.x=1,by.y=col.ID.link)

exp.symb.frame <- data.frame(names(exp.symb2),exp.symb2)
exp.symb.frame <- exp.symb.frame[order(exp.symb.frame[,1]),]

expr_inc_chr_pos_symb <- merge(exp.symb.frame,expr_inc_chr,by.x=1,by.y=col.ID.link)

colnames(expr_inc_chr_pos_symb )[col.ID.link:(col.ID.link+3)] <- c("ID", "Symbol", "CHROMOSOME", "START_POS")
return(expr_inc_chr_pos_symb)
}

