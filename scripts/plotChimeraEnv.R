library(ggplot2)

d=numeric(0)
for (inFile in strsplit(snakemake@input, " ", fixed=T)[[1]]) {
    sample=strsplit(strsplit(inFile, "/")[[1]][2], ".", fixed=T)[[1]][1]
    n=read.table(inFile, as.is=c(2))
    colnames(n) = c("score", "query", "parentA", "parentB", "topParent", "idQM", "idQA", "idQB", "idAB", "idQT", "LY", "LN", "LA", "RY", "RN","RA", "div", "YN")
    n$YN=factor(n$YN, levels=c("Y", "?", "N"))
    n$sample=sample
    n$size = as.numeric(sub(";", "", sapply(strsplit(n$query, "="), function(x) x[4])))
    d = rbind(d, n)
}

cbPalette <- c("red", "darkgrey", "blue")

totals = aggregate(size~sample, d, sum)

a=aggregate(size~sample+YN, d, sum)
colnames(a) = c("sample", "YN", "count")
a$proportion = NA
for (i in 1:length(a$count)) {
    a$proportion[i] = a$count[i] / totals[totals$sample==a$sample[i],]$size
}

write.table(a, snakemake@output[["tab"]], sep="\t", row.names=F)

ggplot(a) + geom_bar(aes(sample, proportion, fill=YN), stat="identity") + scale_y_continuous(labels = scales::percent) + labs(x="", y="proportion of reads", title="Chimera classification in different samples") + scale_fill_manual(values=cbPalette, name=NULL, breaks=c("Y", "N", "?"), labels=c("Chimeric", "Non-chimeric", "Unclear")) + theme_bw()
ggsave(snakemake@output[["relative"]], width=16, height=10)
