library(ggplot2)
d=read.table(snankemake@input[[1]], header=F)
colnames(d) = c("sample", "read", "length", "qual")
m=lm(qual~length+sample, d)

ggplot(d) + geom_point(aes(length, qual, color=sample), alpha=0.5) + geom_smooth(aes(length, qual), method = "lm", linetype = 2)
ggsave(snakemake@output[[1]], width=16, height=10)
