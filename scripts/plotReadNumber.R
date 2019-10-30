library(ggplot2)
library(reshape2)

collectData <- function (set, setname, a_total, a_rel) {
    set$prop=NA
    for (tSample in unique(set$sample)) {
        for (stage in c("lenFilter", "qualFilter", "winQualFilter", "primerFilter")) {
            set[set$sample==tSample & set$stage==stage,]$prop = set[set$sample==tSample & set$stage==stage,]$number/set[set$sample==tSample & set$stage=="raw",]$number
        }
        set[set$sample==tSample & set$stage=="raw",]$prop = 1.0
    }
    new = aggregate(prop~stage, set, mean)
    colnames(new) = c("stage", "mean")
    new$sd=aggregate(prop~stage, set, sd)$prop
    new$group=setname
    a_rel=rbind(a_rel, new)
    new = aggregate(number~stage, set, mean)
    colnames(new) = c("stage", "mean")
    new$sd=aggregate(number~stage, set, sd)$number
    new$group=setname
    a_total=rbind(a_total, new)
    return(list(a_total, a_rel))
}

m=read.table(snakemake@input[[2]], stringsAsFactors=F)
colnames(m) = c("sample", "sampleName", "SRA", "fwdBcId", "fwBcSeq", "revBcId", "reBcSeq", "group")
iso = m[m$group=="isolate",]$sample
env = m[m$group=="env",]$sample

d=read.table(snakemake@input[[1]])
colnames(d) = c("stage", "sample", "sampleName", "number")
d$stage=factor(d$stage, levels=c("raw", "lenFilter", "qualFilter", "winQualFilter", "primerFilter"))
d$group=NA
d[d$sample %in% iso,]$group="isolate"
d[d$sample %in% env,]$group="env"
at=numeric(0)
ar=numeric(0)
env=subset(d, group=="env")
rv=collectData(env, "environmental samples", at, ar)
at=rv[[1]]
ar=rv[[2]]
iso=subset(d, group=="isolate")
rv=collectData(iso, "isolate samples", at, ar)
at=rv[[1]]
ar=rv[[2]]
#ggplot(subset(ar, stage!="raw"), aes(x=stage, y=mean, fill=group)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + labs(y="proportion of reads remaining")
ggplot(subset(ar, stage!="raw"), aes(x=stage, y=mean, fill=group)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + labs(y="proportion of reads remaining") + facet_grid(group~.)
ggsave(snakemake@output[["filtering"]], width=16, height=10)

ggplot(d[d$stage=="raw",], aes(group, number, color=group)) + geom_point(position=position_jitterdodge()) + labs(x="sample group", y="number of reads") + theme(legend.position = "none") + scale_y_log10() + scale_x_discrete(labels=c("environmental samples", "isolates"))
ggsave(snakemake@output[["groups"]], width=16, height=10)

#outtab=dcast(d[!is.na(d$group),], group+sampleName~stage, fun.aggregate=sum, value.var="number")
#outtab=outtab[order(outtab$group, outtab$sampleName),]


