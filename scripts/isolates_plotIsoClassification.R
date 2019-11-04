library(ggplot2)

d=read.table(snakemake@input[[1]], header=T, sep="\t")

d$ssuTax[d$ssuTax=="unknown"] = NA
d$lsuTax[d$lsuTax=="unknown"] = NA
d$itsTax[d$itsTax=="unknown"] = NA

p=ggplot(d) + geom_bar(aes(sampleId, weight=readCount, fill=ssuTax), color="black") + coord_flip() + labs(x="read count", y="sample", fill="SSU classification") + theme_bw()
ggsave(snakemake@output[["ssu"]], p, width=16, height=10)

p=ggplot(d) + geom_bar(aes(sampleId, weight=readCount, fill=itsTax), color="black") + coord_flip() + labs(x="read count", y="sample", fill="ITS classification") + theme_bw()
ggsave(snakemake@output[["its"]], p, width=16, height=10)

p=ggplot(d) + geom_bar(aes(sampleId, weight=readCount, fill=lsuTax), color="black") + coord_flip() + labs(x="read count", y="sample", fill="LSU classification") + theme_bw()
ggsave(snakemake@output[["lsu"]], p, width=16, height=10)

