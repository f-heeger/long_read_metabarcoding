library(ggplot2)

SIZE = snakemake@params[["subset"]]
d=numeric(0)

#inFiles = strsplit("{input}", " ", fixed=T)[[1]]
inFiles = snakemake@input

sampleName = read.table(inFiles[1], row.names=1, as.is=T)

for (inFile in inFiles[2:length(inFiles)]) {
    sample=sampleName[strsplit(strsplit(inFile, "/")[[1]][2], ".", fixed=T)[[1]][1], 1]
    n=read.table(inFile)
    colnames(n) = c("score", "query", "parentA", "parentB", "topParent", "idQM", "idQA", "idQB", "idAB", "idQT", "LY", "LN", "LA", "RY", "RN","RA", "div", "YN")
    n$YN=factor(n$YN, levels=c("Y", "?", "N"))
    n$sample=sample
    d = rbind(d, n)
}

cbPalette <- c("red", "darkgrey", "blue")

ggplot(d) + geom_bar(aes(sample, fill=YN)) + scale_fill_manual(values=cbPalette)
ggsave(snakemake@output[["total"]], width=16, height=10)

totals = aggregate(query~sample, d, length)

a=aggregate(query~sample+YN, d, length)
colnames(a) = c("sample", "YN", "count")
a$proportion = NA
for (i in 1:length(a$count)) {
    a$proportion[i] = a$count[i] / totals[totals$sample==a$sample[i],]$query
}

write.table(a, snakemake@output[["tab"]], sep="\t", row.names=F)

cy = subset(a, sample %in% c("c13i8", "c15i8", "c18i8", "c25i8", "c30i8", "emPCR"))
input = subset(a, sample %in% c("c18i2", "c18i8", "c18i20"))
input$sample = factor(input$sample, levels=c("c18i2", "c18i8", "c18i20"))

ggplot(cy) + geom_bar(aes(sample, proportion, fill=YN), stat="identity") + scale_y_continuous(labels = scales::percent) + labs(x="", y="proportion of reads", title="Chimera classification with different PCR conditions") + scale_x_discrete(labels=c("c13i8"="13 cycles", "c15i8"="15 cyles", "c18i8"="18 cyles", "c25i8"="25 cyles", "c30i8"="30 cyles", "emPCR"="emulsion PCR")) + scale_fill_manual(values=cbPalette, name=NULL, breaks=c("Y", "N", "?"), labels=c("Chimeric", "Non-chimeric", "Unclear")) + theme_bw()
ggsave(snakemake@output[["relative_cy"]], width=16, height=10)

ggplot(input) + geom_bar(aes(sample, proportion, fill=YN), stat="identity") + scale_fill_manual(values=cbPalette) + scale_y_continuous(labels = scales::percent) + labs(x="amount of input template [ng]", y="proportion of reads", title="Chimera classification with different amount of input template") + scale_x_discrete(labels=c("c18i2" = "2", "c18i8"="8", "c18i20"="20"))
ggsave(snakemake@output[["relative_in"]], width=16, height=10)

subData = numeric(0)

for (cond in c("c13i8", "c15i8", "c18i8", "c25i8", "c30i8", "c18i20", "c18i2", "emPCR")) {
    sub = data.frame("N"=rep(NA, 100), "Y"=rep(NA, 100), "?"=rep(NA, 100))
    for (r in 1:100) {
        sampleTab = subset(d, sample==cond)
        sSet = sampleTab[sample(nrow(sampleTab), SIZE),]
        a=aggregate(query~YN, sSet, FUN=length, drop=F)
        
        for (cls in c("N", "?", "Y")) {
            n=a[a$YN==cls,]
            if (dim(n)[1] == 0) {
                sub[r, cls] = 0
            } else {
                sub[r, cls] = n$query
            }
        }
        
    }
    for (cls in c("N", "?", "Y")) {
        subData = rbind(subData, data.frame("sample"=cond, "cls"=cls, "mean"=mean(sub[,cls]), "sd"=sd(sub[,cls])))
    }
}

cy = subset(subData, sample %in% c("c13i8", "c15i8", "c18i8", "c25i8", "c30i8"))
input = subset(subData, sample %in% c("c18i2", "c18i8", "c18i20"))
input$sample = factor(input$sample, levels=c("c18i2", "c18i8", "c18i20"))

ggplot(cy, aes(x=sample, y=mean, ymax=mean+sd, ymin=mean-sd, color=cls)) + geom_point() + geom_linerange() + labs(x="cycles", y="number of reads", title="Chimera classification for different number of PCR cyles \n(subset to 200 reads per sample)") + scale_x_discrete(labels=c("c13i8"="13", "c15i8"="15", "c18i8"="18", "c25i8"="25", "c30i8"="30"))
ggsave(sankemake@output[["subset_cy"]], width=16, height=10)

ggplot(input, aes(x=sample, y=mean, ymax=mean+sd, ymin=mean-sd, color=cls)) + geom_point() + geom_linerange() + labs(x="amount of input template [ng]", y="number of reads", title="Chimera classification for different amount of input template \n(subset to 200 reads per sample)") + scale_x_discrete(labels=c("c18i2" = "2", "c18i8"="8", "c18i20"="20"))
ggsave(snakemake@output[["subset_in"]], width=16, height=10)
