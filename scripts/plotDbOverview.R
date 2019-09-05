library(ggplot2)

colN = c("id", "len", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

ssu=read.table(snakemake@input[["ssu"]], sep="\t")
colnames(ssu) = colN
p1 = ggplot(ssu, aes(x=len)) + geom_histogram(binwidth=1, color="black") + labs(title="Length distribution of SSU data base", x="sequence length", y="count")
#p2 = ggplot(ssu, aes(x=len)) + geom_density()

its=read.table(snakemake@input[["its"]], sep="\t")
colnames(its) = colN
p3 = ggplot(its, aes(x=len)) + geom_histogram(binwidth=1, color="black") + labs(title="Length distribution of ITS data base", x="sequence length", y="count")
#p4 = ggplot(its, aes(x=len)) + geom_density()

lsu=read.table(snakemake@input[["lsu"]], sep="\t")
colnames(lsu) = colN
p5 = ggplot(lsu, aes(x=len)) + geom_histogram(binwidth=1, color="black") + labs(title="Length distribution of LSU data base", x="sequence length", y="count")
#p6 = ggplot(lsu, aes(x=len)) + geom_density()

pdf(snakemake@output[["length"]], width=16, height=10)
print(p1)
#print(p2)
print(p3)
#print(p4)
print(p5)
#print(p6)
dev.off()

ssu$marker="silva (SSU)"
its$marker="UNITE (ITS)"
lsu$marker="RDP Fungal 28S (LSU)"
d=rbind(ssu, its, lsu)
d$marker=factor(d$marker, levels=c("silva (SSU)", "UNITE (ITS)", "RDP Fungal 28S (LSU)"))
p7 = ggplot(d) + geom_bar(aes(x=marker, fill=domain)) + labs(title="Sequences from different domains in the databases", x="data base")
p8 = ggplot(d) + geom_bar(aes(x=marker, fill=kingdom)) + labs(title="Sequences from different kingdoms in the databases", x="data base")
p9 = ggplot(d[d$kingdom=="Fungi",]) + geom_bar(aes(x=marker, fill=phylum)) + labs(title="Sequences from different fungal kingdoms in the databases", x="data base")
pdf(snakemake@output[["tax"]], width=16, height=10)
print(p7)
print(p8)
print(p9)
dev.off()
