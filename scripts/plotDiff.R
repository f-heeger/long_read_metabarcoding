library(ggplot2)
library(reshape2)

d=read.table(snakemake@input[[1]], sep="\t")
colnames(d) = c("oId", "size", "fungi", "mrk1", "mrk2", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

for (lvl in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
    d[,lvl] = factor(d[,lvl], levels=c("same", "diff", "blank", NA))
}

m=melt(d, id.vars=c("oId", "size", "fungi", "mrk1", "mrk2"))
m$mrk1=factor(m$mrk1, levels=c("ssu", "its", "lsu"))
m$mrk2=factor(m$mrk2, levels=c("ssu", "its", "lsu"))
m$value = factor(m$value, levels=c(NA, "blank", "diff", "same"))

m$mrk1=factor(m$mrk1, levels=c("ssu", "its", "lsu"))
m$mrk2=factor(m$mrk2, levels=c("ssu", "its", "lsu"))
ggplot(subset(m, size>1 & !is.na(value) & variable!="domain")) + geom_bar(aes(x=variable, fill=value)) + facet_grid(mrk1~mrk2) + scale_fill_manual(values=c("grey30", "red", "blue")) + labs(y="OTU count", x="taxonomic level") + theme_bw()

#ggplot(subset(m, size>1 & !is.na(value) & variable!="domain")) + geom_bar(aes(x=variable, fill=value)) + facet_grid(mrk1~mrk2, switch="y") + scale_fill_manual(values=c("grey30", "red", "blue")) + labs(y="OTU count", x="taxonomic level") + scale_y_continuous(position = "right") + theme_bw()

ggsave(snakemake@output[["bars"]], width=16, height=10)


#proportion of comparable (assigned in both markers) OTUs
mrk = levels(m$mrk1)
tax = levels(m$variable)

a=aggregate(oId~mrk1+mrk2+variable+value, m[m$size>1,], length)
a$prop = NA
a$pAss = NA
for (m1 in mrk) {
    for (m2 in setdiff(mrk, c(m1))) {
        for (v1 in tax) {
            total = sum(subset(a, mrk1==m1 & mrk2==m2 & variable==v1)$oId)
            ass = sum(subset(a, mrk1==m1 & mrk2==m2 & variable==v1 & (value=="same" | value=="diff"))$oId)
            for (v2 in unique(subset(a, mrk1==m1 & mrk2==m2 & a$variable==v1)$value)) {
                a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$prop = a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$oId/total
            }
            for (v2 in setdiff(unique(subset(a, mrk1==m1 & mrk2==m2 & a$variable==v1)$value), c("blank"))) {
                a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$pAss = a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$oId/ass
            }
        }
    }
}

#ggplot(a) + geom_line(aes(x=variable, y=pAss, group=value, color=value)) + facet_grid(mrk1~mrk2, switch="y") + scale_color_manual(values=c("grey30", "red", "blue")) + labs(y="proportion of comparable OTUs", x="taxonomic level") + theme_bw()
a$comp=paste(a$mrk1, a$mrk2, sep="-")
ggplot(subset(a, value=="same" & (comp %in% c("ssu-its", "ssu-lsu", "its-lsu")) & variable!="domain")) + geom_point(aes(x=variable, y=pAss, color=comp)) + geom_line(aes(x=variable, y=pAss, group=comp, color=comp), linetype=2) + theme_bw() + labs(y="proportion of comparable OTUs that have the same classification", x="taxonomic level")
ggsave(snakemake@output[["prop"]], width=16, height=10)
