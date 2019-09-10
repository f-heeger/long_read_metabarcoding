d=read.table(snakemake@input[[1]], header=T)

#nodes
nodes=data.frame(id=numeric(0), size=numeric(0), col=numeric(0))
for (type in c("denovo", "ref")) {
    n=aggregate(as.formula(paste("readId~", type, sep="")), d, FUN=length)
    for (i in 1:dim(n)[1]) {
        c="grey"
        if (n[i,1] == "N") {
            c="blue"
        }
        if (n[i,1] == "Y") {
            c="red"
        }
        nodes=rbind(nodes, data.frame(id=paste(type, n[i,1], sep="_"), size=n[i,2], col=c, stringsAsFactors=F))
    }
}

e=aggregate(readId~denovo+ref, d, FUN=length)
edges=data.frame(from=paste("denovo", e$denovo, sep="_"), to=paste("ref", e$ref, sep="_"), weight=e$readId, stringsAsFactors=F)


library(sankey)
cls_sankey <- make_sankey(nodes=nodes, edges = edges)
svg(snakemake@output[[1]], width=16, height=10)
sankey(cls_sankey)
dev.off()
