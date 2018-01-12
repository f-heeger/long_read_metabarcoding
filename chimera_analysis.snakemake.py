import pickle
import random

from snakemake.utils import min_version, R

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

#inFiles = {
#"emPCR" : "../PacBioMetabarcoding2/raw/0009_Forward--0009_Forward.fastq",
#"c13i8": "../PacBioMetabarcoding2/raw/0018_Forward--0018_Forward.fastq",
#"c15i8": "../PacBioMetabarcoding2/raw/0027_Forward--0027_Forward.fastq",
#"c18i20": "../PacBioMetabarcoding2/raw/0034_Forward--0034_Forward.fastq",
#"c18i8": "../PacBioMetabarcoding2/raw/0056_Forward--0056_Forward.fastq",
#"c25i8": "../PacBioMetabarcoding2/raw/0075_Forward--0075_Forward.fastq",
#"c18i2": "../PacBioMetabarcoding2/raw/0095_Forward--0095_Forward.fastq",
#"c30i8" : "../PacBioMetabarcoding3b/raw/Lib4-0018.fastq"
#}

#configfile: "config.json"

#shell.prefix("sleep 3; ") #work around to desl with "too quck" rule execution and slow samba

#rule all:
#    input: "chimeraCyclesRelativeBarplot.svg", "chimera_comp_sankey.svg"

#def getInfilePath(wildcards):
#    return inFiles[wildcards.sample]

#rule lengthFilter:
#    input: getInfilePath
#    output: right="lenFilter/{sample}_rightLen.fastq", long="lenFilter/{sample}_tooLong.fastq", short="lenFilter/{sample}_tooShort.fastq"
#    log: "logs/{sample}_lenFilter.log"
#    params: maxLen = 6500, minLen = 3000
#    run:
#        nrLong = 0
#        nrShort = 0
#        with open(output.right, "w") as out, open(output.long, "w") as tooLong, open(output.short, "w") as tooShort:
#            for read in SeqIO.parse(open(input[0]), "fastq"):
#                if len(read) <= params.minLen:
#                    tooShort.write(read.format("fastq"))
#                    nrShort += 1
#                elif len(read) <= params.maxLen:
#                    out.write(read.format("fastq"))
#                else:
#                    tooLong.write(read.format("fastq"))
#                    nrLong += 1
#        with open(log[0], "w") as logFile:
#            logFile.write("%i reads were removed because they were longer than %i\n" % (nrLong, params.maxLen))
#            logFile.write("%i reads were removed because they were shorter than %i\n" % (nrLong, params.maxLen))

#rule qualityFilter:
#    input: fastq="lenFilter/{sample}_rightLen.fastq"
#    output: good="qualFilter/{sample}_goodQual.fastq", bad="qualFilter/{sample}_badQual.fastq", info="qualFilter/{sample}_qualInfo.tsv"
#    params: minQual=0.996
#    log: "logs/{sample}_qualityFilter.log"
#    run:
#        removed = 0
#        with open(output.good, "w") as out, open(output.bad, "w") as trash, open(output.info, "w") as info:
#            for read in SeqIO.parse(open(input.fastq), "fastq"):
#                try:
#                    tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"]]) / len(read)
#                except Exception:
#                    print(read.id)
#                    raise
#                info.write("%s\t%f\n" % (read.id, 1-tError))
#                if (1-tError) < params.minQual:
#                    removed += 1
#                    trash.write(read.format("fastq"))
#                else:
#                    out.write(read.format("fastq"))
#        open(log[0], "w").write("%s: %i reads removed because quality < %f\n" % (wildcards.sample, removed, params.minQual))

#rule windowQualFilter:
#    input: fastq="qualFilter/{sample}_goodQual.fastq"
#    output: good="windowQualFilter/{sample}_goodQual.fastq"
#    params: winSize=8, minQual=0.9
#    log: "logs/{sample}_winQualityFilter.log"
#    run:
#        removed = 0
#        with open(output.good, "w") as out:
#            for read in SeqIO.parse(open(input.fastq), "fastq"):
#                tRemoved=False
#                for i in range(len(read)-params.winSize):
#                    tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"][i:i+params.winSize] ]) / params.winSize
#                    if (1-tError) < params.minQual:
#                        removed += 1
#                        tRemoved = True
#                        break
#                if not tRemoved:
#                    out.write(read.format("fastq"))
#            open(log[0], "w").write("%s: %i reads removed because in a window of size %i quality droped below %f\n" % (wildcards.sample, removed, params.winSize, params.minQual))


#rule init_filterPrimer:
#    input: "windowQualFilter/{sample}_goodQual.fastq"
#    output: fastq="primers/{sample}_primer.fastq"
#    log: "logs/{sample}_primer.log"
#    threads: 3
#    params: minovl=10
#    run:
#        shell("%(cutadapt)s -g forward=%(fwd_primer)s -g reverse=%(rev_primer)s -O {params.minovl} --trimmed-only -o primers/{wildcards.sample}_{{name}}.fastq {input} &> {log}" % config)
#        config["fwd_primer_rc"] = str(Seq(config["fwd_primer"], IUPACAmbiguousDNA()).reverse_complement())
#        config["rev_primer_rc"] = str(Seq(config["rev_primer"], IUPACAmbiguousDNA()).reverse_complement())
#        
#        shell("%(cutadapt)s -a forward=%(fwd_primer_rc)s -O {params.minovl} --trimmed-only -o primers/{wildcards.sample}_reverse_forward.fastq primers/{wildcards.sample}_reverse.fastq &>> {log}" % config)
#        shell("%(cutadapt)s -a reverse=%(rev_primer_rc)s -O {params.minovl} --trimmed-only -o primers/{wildcards.sample}_forward_reverse.fastq primers/{wildcards.sample}_forward.fastq &>> {log}" % config)
#        with open(output[0], "w") as out:
#            for rec in SeqIO.parse(open("primers/%s_reverse_forward.fastq" % wildcards.sample), 
#                                   "fastq"):
#                newRec = rec.reverse_complement()
#                newRec.id = "%s/rc" % rec.id
#                newRec.description = rec.description.split(" ", 1)[1]
#                out.write(newRec.format("fastq"))
#        shell("cat primers/{wildcards.sample}_forward_reverse.fastq >> {output}")

#rule fastq2fasta:
#    input: "primers/{sample}_primer.fastq"
#    output: "primers/{sample}_primer.fasta"
#    run:
#        with open(output[0], "w") as out:
#            for read in SeqIO.parse(open(input[0]), "fastq"):
#                out.write(read.format("fasta"))

#rule prepPrecluster:
#    input: fasta="primers/{sample}_primer.fasta", qual="qualFilter/{sample}_qualInfo.tsv"
#    output: "preclusters/{sample}_cluInput.fasta"
#    run:
#        qual = {}
#        for line in open(input.qual):
#            rId, tQual = line.strip().split("\t")
#            qual[rId] = tQual
#        readList = []
#        for read in SeqIO.parse(open(input.fasta), "fasta"):
#            if read.id.endswith("rc"):
#                readId = read.id.rsplit("/", 1)[0]
#            else:
#                readId = read.id
#            readList.append((read, qual[readId]))
#        with open(output[0], "w") as out:
#            for read, qual in sorted(readList, key=lambda x: x[1], reverse=True):
#                read.description=str(qual)
#                out.write(read.format("fasta"))

#rule preCluster:
#    input: "preclusters/{sample}_cluInput.fasta"
#    output: uc="preclusters/{sample}.uc.txt", cons="consensus/{sample}_consensus.fasta"
#    log: "logs/{sample}_pre-cluster.log"
#    threads: 6
#    shell:
#        "%(vsearch)s --usersort --cluster_smallmem {input} --relabel precluster --sizeout --iddef 0 --id 0.99 --minsl 0.9 --consout {output.cons} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

#rule preClusterReads:
#    input: clsInfo="preclusters/{sample}.uc.txt"
#    output: info="preclusters/{sample}_cluster_reads.pic"
#    run:
#        cluster2read = {}
#        for line in open(input.clsInfo):
#            if line[0] == "C":
#                pass
#            elif line[0] in "SH":
#                arr = line.strip().split("\t")
#                seq, clu = arr[-2:]
#                if clu=="*":
#                    cluster2read[seq] = [seq]
#                else:
#                    cluster2read[clu].append(seq)
#            else:
#                raise ValueError("Unknown record type: %s" % arr[0])
#        pickle.dump(cluster2read, open(output.info, "wb"))

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

rule plotChimera:
    input: expand("refChimera/{sample}.chimeraReport.tsv", sample=mockSamples)
    output: total="chimeraTotalBarplot.svg", tab="chimeraTable.tsv", relative_cy="chimeraCyclesRelativeBarplot.svg", relative_in="chimeraInputRelativeBarplot.svg", subset_cy="chimeraCyclesSubset.svg", subset_in="chimeraInputSubset.svg"
    params: subset=200
    run:
        R("""
        library(ggplot2)
        SIZE = {params.subset}
        d=numeric(0)

        for (inFile in strsplit("{input}", " ", fixed=T)[[1]]) {{
            sample=strsplit(strsplit(inFile, "/")[[1]][2], ".", fixed=T)[[1]][1]
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
    input: denovo="denovoChimera/{sample}.chimeraReport.tsv", ref="refChimera{sample}.chimeraReport.tsv", cluster2read="preclusters/{sample}_cluInfo.tsv"
    output: "refChimera{sample}_comp.tsv"
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

#rule compareChimera:
#    input: denovo="denovoChimera/{sample}.chimeraReport.tsv", ref="chimera/{sample}.chimeraReport.tsv", cluster2read="preclusters/{sample}_cluster_reads.pic"
#    output: "chimera/{sample}_comp.tsv"
#    run:
#        clu2read = pickle.load(open(input.cluster2read, "rb"))
#        ref={}
#        for line in open(input.ref):
#            arr = line.strip().split("\t")
#            scr=arr[0]
#            readId = arr[1]
#            cls = arr[-1]
#            ref[readId] = (cls, scr)
#        with open(output[0], "w") as out:
#            out.write("readId\tclusterNr\tclusterId\tdenovoScr\tdenovo\trefScr\tref\n")
#            for l,line in enumerate(open(input.denovo)):
#                arr = line.strip().split("\t")
#                scr=arr[0]
#                cluId = arr[1].split("=", 1)[1].split(";", 1)[0]
#                cls = arr[-1]
#                for rId in clu2read[cluId]:
#                    out.write("%s\tcluster%i\t%s\t%s\t%s\t%s\t%s\n" % (rId, l, cluId, scr, cls, ref[rId][1], ref[rId][0]))

rule plotChimeraComp:
    input: "refChimerac30i8_comp.tsv"
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

