library(ggplot2)
library(reshape2)

d=read.table(snakemake@input[[1]], header=T)
d$sample = factor(d$sample, levels=d$sample)

abs = d
abs$totalRemoved = rowSums(abs[,3:dim(abs)[2]])
abs$remaining = abs$raw - abs$totalRemoved

for (i in 3:dim(d)[2]) {
    d[,i] = d[,i]/d[,2]
}


m = melt(d, id.vars=c("sample"))

#labs=data.frame(text=abs$remaining, x=1:dim(abs)[1], y=rowSums(d[,3:dim(d)[2]])+0.05)
labs = data.frame(text=abs$remaining, x=1:dim(abs)[1], y=rep(-0.05,dim(abs)[1]))

ggplot(m[m$variable!="raw",]) + geom_bar(aes(sample, fill = variable, weight=value)) + ylim(-0.05,1) + annotate("text", label=labs$text, x=labs$x, y=labs$y) + coord_flip() + labs(fill="removal reason", y="proportion of reads", title="Proportion of reads removed by filtering steps\n (number on the left are total remaining reads)") + scale_fill_brewer(type="qual", palette=6, labels=c("too long", "too short", "low mean quality", "low window quality", "no primer found", "chimera", "no ITS found"))

ggsave(snakemake@output[[1]])
