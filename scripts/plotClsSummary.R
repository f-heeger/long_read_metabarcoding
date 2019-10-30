library(reshape2)
library(ggplot2)

d=read.table(snakemake@input[["all"]], sep="\t")
colnames(d) = c("oId", "size", "marker", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "depth")

d$marker=factor(d$marker, levels=c("ssu", "its", "lsu"))
d$oId=factor(d$oId, levels=unique(d[order(d$size, decreasing=T),]$oId))
ggplot(d[1:600,], aes(x=oId,fill=phylum, weight=depth)) + geom_bar() + facet_grid(marker~.) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
ggsave(snakemake@output[["depth"]], width=16, height=10)

fo = unique(d[d$kingdom=="Fungi",]$oId)
f=subset(d, d$oId %in% fo)

f$oId=factor(f$oId, levels=unique(f[order(f$marker, f$phylum, -f$depth),]$oId))
 ggplot(f[1:600,], aes(x=oId,fill=phylum, weight=depth)) + geom_bar() + facet_grid(marker~.) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6)) + scale_y_continuous(breaks=seq(1, 7, 1), labels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
ggsave(snakemake@output[["depthFungi"]], width=16, height=10)

#sum(table(subset(f[1:600,], marker=="ssu" & size>1)$depth)[c("2","3","4","5","6","7")])
#sum(table(subset(f[1:600,], marker=="lsu" & size>1)$depth)[c("2","3","4","5","6","7")])
#sum(table(subset(f[1:600,], marker=="its" & size>1)$depth)[c("2","3","4","5","6","7")])

#d$marker=factor(d$marker, levels=c("lsu", "ssu", "its"))
#ggplot(d[1:600,], aes(x=oId, y=marker, fill=phylum)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
#ggsave(snakemake@output[["block"]], width=16, height=10)
