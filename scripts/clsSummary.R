d=read.table(snakemake@input[[1]], sep="\t")
colnames(d) = c("oId", "size", "mrk", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "depth")
a=aggregate(depth ~ oId, subset(d, size>1), max)
ranks=c("kingdom", "phylum", "class", "order", "family", "genus", "species")
s = data.frame(rank=numeric(0), nubmer=numeric(0))
for (i in 1:7) {
    s = rbind(s, data.frame(rank=ranks[i], number=sum(a$depth>=i)))
}
write.table(s, snakemake@output[[1]], sep="\t", row.names=F)
