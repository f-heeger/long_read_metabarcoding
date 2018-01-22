import random

rule removeChimeraRef:
    input: seqs="primers/{sample}_primer.fasta", ref="isolateSeqs.fasta"
    output: fasta="refChimera/{sample}.nochimera.fasta", tsv="refChimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_refChimera.log"
    shell:
        "%(vsearch)s --uchime_ref {input.seqs} --db {input.ref} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule subset:
    input: seqs="primers/{sample}_primer.fasta", tsv="refChimera/{sample}.chimeraReport.tsv"
    output: "manualCheck/{sample}_subset.fasta"
    params: N=100
    run:
        random.seed("42")
        chimera = {}
        for line in open(input.tsv):
            arr=line.strip().rsplit("\t")
            chimera[arr[1]] = arr[-1]
        reads = list(SeqIO.parse(open(input.seqs), "fasta"))
        subset = random.sample(reads, params.N)
        with open(output[0], "w") as out:
            for r, read in enumerate(subset):
                read.description = chimera[read.id]
                read.id = "%i_%s" % (r, read.id)
                out.write(read.format("fasta"))

rule sampleNamesTab:
    output: "mockSamples.tsv"
    run:
        with open(output[0], "w") as out:
            for sId, sample in sampleInfo["sampleName"].items():
                if sampleInfo["sampleType"][sId] == "mock":
                    out.write("%s\t%s\n" % (sId, sample))

rule plotChimera:
    input: "mockSamples.tsv", expand("refChimera/{sample}.chimeraReport.tsv", sample=mockSamples)
    output: total="chimeraTotalBarplot.svg", tab="chimeraTable.tsv", relative_cy="chimeraCyclesRelativeBarplot.svg", relative_in="chimeraInputRelativeBarplot.svg", subset_cy="chimeraCyclesSubset.svg", subset_in="chimeraInputSubset.svg"
    params: subset=200
    run:
        R("""
        library(ggplot2)
        SIZE = {params.subset}
        d=numeric(0)

        inFiles = strsplit("{input}", " ", fixed=T)[[1]]
        
        sampleName = read.table(inFiles[1], row.names=1, as.is=T)
        
        for (inFile in inFiles[2:length(inFiles)]) {{
            sample=sampleName[strsplit(strsplit(inFile, "/")[[1]][2], ".", fixed=T)[[1]][1], 1]
            n=read.table(inFile)
            colnames(n) = c("score", "query", "parentA", "parentB", "topParent", "idQM", "idQA", "idQB", "idAB", "idQT", "LY", "LN", "LA", "RY", "RN","RA", "div", "YN")
            n$YN=factor(n$YN, levels=c("Y", "?", "N"))
            n$sample=sample
            d = rbind(d, n)
        }}

        cbPalette <- c("red", "darkgrey", "blue")

        ggplot(d) + geom_bar(aes(sample, fill=YN)) + scale_fill_manual(values=cbPalette)
        ggsave("{output.total}", width=16, height=10)

        totals = aggregate(query~sample, d, length)

        a=aggregate(query~sample+YN, d, length)
        colnames(a) = c("sample", "YN", "count")
        a$proportion = NA
        for (i in 1:length(a$count)) {{
            a$proportion[i] = a$count[i] / totals[totals$sample==a$sample[i],]$query
        }}

        write.table(a, "{output.tab}", sep="\t", row.names=F)

        cy = subset(a, sample %in% c("c13i8", "c15i8", "c18i8", "c25i8", "c30i8", "emPCR"))
        input = subset(a, sample %in% c("c18i2", "c18i8", "c18i20"))
        input$sample = factor(input$sample, levels=c("c18i2", "c18i8", "c18i20"))

        ggplot(cy) + geom_bar(aes(sample, proportion, fill=YN), stat="identity") + scale_y_continuous(labels = scales::percent) + labs(x="", y="proportion of reads", title="Chimera classification with different PCR conditions") + scale_x_discrete(labels=c("c13i8"="13 cycles", "c15i8"="15 cyles", "c18i8"="18 cyles", "c25i8"="25 cyles", "c30i8"="30 cyles", "emPCR"="emulsion PCR")) + scale_fill_manual(values=cbPalette, name=NULL, breaks=c("Y", "N", "?"), labels=c("Chimeric", "Non-chimeric", "Unclear")) + theme_bw()
        ggsave("{output.relative_cy}", width=16, height=10)

        ggplot(input) + geom_bar(aes(sample, proportion, fill=YN), stat="identity") + scale_fill_manual(values=cbPalette) + scale_y_continuous(labels = scales::percent) + labs(x="amount of input template [ng]", y="proportion of reads", title="Chimera classification with different amount of input template") + scale_x_discrete(labels=c("c18i2" = "2", "c18i8"="8", "c18i20"="20"))
        ggsave("{output.relative_in}", width=16, height=10)

        subData = numeric(0)

        for (cond in c("c13i8", "c15i8", "c18i8", "c25i8", "c30i8", "c18i20", "c18i2", "emPCR")) {{
            sub = data.frame("N"=rep(NA, 100), "Y"=rep(NA, 100), "?"=rep(NA, 100))
            for (r in 1:100) {{
                sampleTab = subset(d, sample==cond)
                sSet = sampleTab[sample(nrow(sampleTab), SIZE),]
                a=aggregate(query~YN, sSet, FUN=length, drop=F)
                
                for (cls in c("N", "?", "Y")) {{
                    n=a[a$YN==cls,]
                    if (dim(n)[1] == 0) {{
                        sub[r, cls] = 0
                    }} else {{
                        sub[r, cls] = n$query
                    }}
                }}
                
            }}
            for (cls in c("N", "?", "Y")) {{
                subData = rbind(subData, data.frame("sample"=cond, "cls"=cls, "mean"=mean(sub[,cls]), "sd"=sd(sub[,cls])))
            }}
        }}

        cy = subset(subData, sample %in% c("c13i8", "c15i8", "c18i8", "c25i8", "c30i8"))
        input = subset(subData, sample %in% c("c18i2", "c18i8", "c18i20"))
        input$sample = factor(input$sample, levels=c("c18i2", "c18i8", "c18i20"))

        ggplot(cy, aes(x=sample, y=mean, ymax=mean+sd, ymin=mean-sd, color=cls)) + geom_point() + geom_linerange() + labs(x="cycles", y="number of reads", title="Chimera classification for different number of PCR cyles \n(subset to 200 reads per sample)") + scale_x_discrete(labels=c("c13i8"="13", "c15i8"="15", "c18i8"="18", "c25i8"="25", "c30i8"="30"))
        ggsave("{output.subset_cy}", width=16, height=10)

        ggplot(input, aes(x=sample, y=mean, ymax=mean+sd, ymin=mean-sd, color=cls)) + geom_point() + geom_linerange() + labs(x="amount of input template [ng]", y="number of reads", title="Chimera classification for different amount of input template \n(subset to 200 reads per sample)") + scale_x_discrete(labels=c("c18i2" = "2", "c18i8"="8", "c18i20"="20"))
        ggsave("{output.subset_in}", width=16, height=10)
        """)

rule removeChimeraDenovo:
    input: seqs="consensus/{sample}_consensus.fasta"
    output: fasta="denovoChimera/{sample}.nochimera.fasta", tsv="denovoChimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_denovoChimera.log"
    shell:
        "%(vsearch)s --uchime_denovo {input.seqs} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule compareChimera:
    input: denovo="denovoChimera/{sample}.chimeraReport.tsv", ref="refChimera/{sample}.chimeraReport.tsv", cluster2read="preclusters/{sample}_cluInfo.tsv"
    output: "refChimera_{sample}_comp.tsv"
    run:
        clu2read = {}
        for line in open(input.cluster2read):
            clu, read = line.strip().split("\t")
            try:
                clu2read[clu].append(read)
            except KeyError:
                clu2read[clu] = [read]
        ref={}
        for line in open(input.ref):
            arr = line.strip().split("\t")
            scr=arr[0]
            readId = arr[1]
            cls = arr[-1]
            ref[readId] = (cls, scr)
        with open(output[0], "w") as out:
            out.write("readId\tclusterNr\tclusterId\tdenovoScr\tdenovo\trefScr\tref\n")
            for l,line in enumerate(open(input.denovo)):
                arr = line.strip().split("\t")
                scr=arr[0]
                cluId = arr[1].split("=", 1)[1].split(";", 1)[0]
                cls = arr[-1]
                for rId in clu2read[cluId]:
                    out.write("%s\tcluster%i\t%s\t%s\t%s\t%s\t%s\n" % (rId, l, cluId, scr, cls, ref[rId][1], ref[rId][0]))


rule plotChimeraComp:
    input: "refChimera_Lib4-0018_comp.tsv"
    output: "chimera_comp_sankey.svg"
    run:
        R("""
        d=read.table("{input}", header=T)

        #nodes
        nodes=data.frame(id=numeric(0), size=numeric(0), col=numeric(0))
        for (type in c("denovo", "ref")) {{
            n=aggregate(as.formula(paste("readId~", type, sep="")), d, FUN=length)
            for (i in 1:dim(n)[1]) {{
                c="grey"
                if (n[i,1] == "N") {{
                    c="blue"
                }}
                if (n[i,1] == "Y") {{
                    c="red"
                }}
                nodes=rbind(nodes, data.frame(id=paste(type, n[i,1], sep="_"), size=n[i,2], col=c, stringsAsFactors=F))
            }}
        }}

        e=aggregate(readId~denovo+ref, d, FUN=length)
        edges=data.frame(from=paste("denovo", e$denovo, sep="_"), to=paste("ref", e$ref, sep="_"), weight=e$readId, stringsAsFactors=F)


        library(sankey)
        cls_sankey <- make_sankey(nodes=nodes, edges = edges)
        svg("{output}", width=16, height=10)
        sankey(cls_sankey)
        dev.off()

        #library(ggplot2)
        #ggplot(d) + geom_bar(aes(x=clusterNr, fill=ref)) + facet_grid(denovo~.) + scale_y_log10()
        """)

