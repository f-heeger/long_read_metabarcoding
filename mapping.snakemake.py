from snakemake.utils import min_version, R
from Bio import SeqIO

min_version("3.5.4")

shell.prefix("sleep 10; ")

configfile: "config.json"

samples = {"Lib4-0018": "mock",
           "Lib1-0009": "CA1",
           "Lib1-0027": "PC1",
           "Lib1-0056": "Csp1",
           "Lib1-0095": "CL1",
           "Lib2-0009": "UM1",
           "Lib2-0027": "PC2",
           "Lib2-0056": "ME1",
           "Lib3-0009": "SC1",
           "Lib3-0027": "LS1",
           "Lib3-0056": "Csp2",
           "Lib3-0095": "DT1",
           "Lib4-0009": "MR1",
           "Lib4-0027": "EV1",
           "Lib4-0056": "ME2",
           "Lib5-0009": "CA2",
           "Lib5-0027": "EV2",
           "Lib5-0095": "DT2",
           "Lib6-0009": "UM2",
           #"Lib6-0027": "IF1",
           "Lib6-0095": "CL2",
           "Lib7-0009": "SC2",
           "Lib7-0027": "LS2",
           "Lib7-0056": "PB1",
           "Lib7-0095": "Psp1",
           "Lib8-0009": "MR2",
#           "Lib8-0027": "IF2",
           "Lib8-0056": "PB2",
           "Lib8-0095": "Psp2",
           "Lib0-0009": "CHY1",
           "Lib0-0075": "CR",
           "Lib0-0056": "TR"
           }

isolates = {"CA": ["Lib1-0009", "Lib5-0009"],
            "SC": ["Lib3-0009", "Lib7-0009"],
            "CL": ["Lib1-0095", "Lib6-0095"],
            "EV": ["Lib4-0027", "Lib5-0027"],
            "PB": ["Lib7-0056", "Lib8-0056"],
            "CHY1": ["Lib0-0009"],
            "TR": ["Lib0-0056"], 
            "PC": ["Lib1-0027", "Lib2-0027"],
            "Csp": ["Lib1-0056", "Lib3-0056"],
            "Psp": ["Lib7-0095", "Lib8-0095"],
            "CR": ["Lib0-0075"],
            "ME": ["Lib2-0056", "Lib4-0056"],
            "UM": ["Lib2-0009", "Lib2-0009"],
            "LS": ["Lib3-0027", "Lib7-0027"],
            "DT": ["Lib3-0095", "Lib5-0095"],
            #"IF": ["Lib6-0027", "Lib8-0027"],
            "MR": ["Lib4-0009", "Lib8-0009"]
            }

ref = ["LSU", "SSU", "ITS"]

rule mapping:
    input: reads="primers/{sample}_primer.fasta", ref="../PacBioMetabarcoding2/references/all_{ref}.fasta"
    output: m5="mapping/mapping/{sample}_vs_{ref}.m5"
    threads: 3
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule getCls:
    input: "mapping/mapping/{sample}_vs_{ref}.m5"
    output: "mapping/matches/match_{sample}_{ref}.tsv"
    run:
        data={}

        for line in open(input[0]):
            qName, qLength, qStart, qEnd, qStrand, tName, tLength, tStart, tEnd, tStrand, score, numMatch, numMismatch, numIns, numDel, mapQV, qAlignedSeq, matchPattern, tAlignedSeq = [x for x in line.split(" ") if len(x)>0]
            tcov = (int(tEnd)-int(tStart))/float(tLength)
            if tcov < 0.9:
                continue
            ident = float(numMatch)/(int(qEnd)-int(qStart))
            spec, marker = tName.split("/")[0].rsplit("_", 1)
            assert marker == wildcards.ref
            mData = (spec, ident, float(score))
            try:
                data[qName].append(mData)
            except KeyError:
                data[qName] = [mData]
            
        with open(output[0], "w") as out:
            for qName, entryList in data.items():
                if len(entryList) == 1:
                    out.write(qName)
                    out.write("\t%s\t%s\t%s\n" % entryList[0])
                    continue
                entryList.sort(key=lambda x: x[1])
                best = entryList[0]
                for entry in entryList[1:]:
                    iDiff = best[1] - entry[1]
                    if iDiff > 0.01:
                        #this entry is 1% points less similar than the best: stop here
                        break
                    else:
                        out.write(qName)
                        out.write("\t%s\t%s\t%s\n" % entry)
                
                
rule compareMappingCls:
    input: reads="primers/{sample}_primer.fasta", its="mapping/matches/match_{sample}_ITS.tsv", ssu="mapping/matches/match_{sample}_SSU.tsv", lsu="mapping/matches/match_{sample}_LSU.tsv"
    output: "mapping/{sample}_cls.tsv"
    run:
        its = {}
        for line in open(input.its):
            read, cls, ident, scr = line.strip().split("\t")
            try:
                its[read].append(cls)
            except KeyError:
                its[read] = [cls]
        ssu = {}
        for line in open(input.ssu):
            read, cls, ident, scr = line.strip().split("\t")
            try:
                ssu[read].append(cls)
            except KeyError:
                ssu[read] = [cls]
        lsu = {}
        for line in open(input.lsu):
            read, cls, ident, scr = line.strip().split("\t")
            try:
                lsu[read].append(cls)
            except KeyError:
                lsu[read] = [cls]
        
        unclear=0
        unknown=0
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input.reads), "fasta"):
                comb = set(its.get(rec.id,[])) & set(ssu.get(rec.id,[])) & set(lsu.get(rec.id,[]))
                if len(comb) == 1:
                    cls = comb.pop()
                elif len(comb) == 0:
                    cls = "unknown"
                    unknown+=1
                else:
                    cls = "|".join(comb)
                    unclear+=1
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (samples[wildcards.sample],
                                                        rec.id, 
                                                        "|".join(set(ssu.get(rec.id,["unknown"]))), 
                                                        "|".join(set(its.get(rec.id,["unknown"]))), 
                                                        "|".join(set(lsu.get(rec.id,["unknown"]))),
                                                        cls))
        print("%s: %i unknown, %i unclear" % (wildcards.sample, unknown, unclear))
        
rule concatCls:
    input: expand("mapping/{sample}_cls.tsv", sample=samples)
    output: "mapping/combinedCls.tsv"
    shell:
        "cat {input} > {output}"
        
rule plot:
    input: "mapping/combinedCls.tsv"
    output: bar="mapping/readPerSampleBarplot.pdf", heat="mapping/readPerSampleHeatmap.pdf", tab="mapping/readClassification.tsv"
    run:
        R("""
        library(ggplot2)

        d=read.table("{input}", sep="\t")
        colnames(d) = c("sample", "read", "ssu_cls", "its_cls", "lsu_cls", "cls")

        d$cls[d$cls=="unknown"] = NA
        ggplot(d) + geom_bar(aes(sample, fill=cls)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output.bar}")
        a=aggregate(. ~ sample + cls, d, length)
        ggplot(a) + geom_raster(aes(sample, cls, fill=log(read))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output.heat}")
        a=aggregate(read~cls,d,FUN=length)
        colnames(a) = c("cls", "readNumber")
        write.table(a, "{output.tab}", quote=F, sep="\t", row.names=F)
        """)

rule seqByCls:
    input: reads="primers/{sample}_primer.fasta", cls="mapping/combinedCls.tsv"
    output: dynamic("mapping/byCls/{sample}-{cls}.fasta")
    run:
        cls = {}
        outStreams = {}
        try:
            for line in open(input.cls):
                smpl, readId, ssu, its, lsu, comb = line.strip().split("\t")
                cls[readId] = comb
                if comb not in outStreams:
                    outStreams[comb] = open("mapping/byCls/%s-%s.fasta" % (wildcards.sample, comb), "w")
            for rec in SeqIO.parse(open(input.reads), "fasta"):
                outStreams[cls[rec.id]].write(rec.format("fasta"))
        finally:
            for out in outStreams.values():
                out.close()

rule byClsAlignemnet:
    input: "mapping/byCls/{sample}-{cls}.fasta"
    output: "mapping/alnByCls/{sample}-{cls}.aln.fasta"
    log: "mapping/logs/byClsAln_{sample}-{cls}.log"
    threads: 4
    shell:
        "%(mafft)s --ep 1 --thread {threads} {input} > {output} 2> {log}" % config

        
rule clusterMapping:
    input: reads="clusters/{sample}_cluster.fasta", ref="../PacBioMetabarcoding2/references/all_{ref}.fasta"
    output: m5="mapping/cluMapping/{sample}_vs_{ref}.m5"
    threads: 3
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule getCluCls:
    input: "mapping/cluMapping/{sample}_vs_{ref}.m5"
    output: "mapping/cluMatches/match_{sample}_{ref}.tsv"
    run:
        data={}

        for line in open(input[0]):
            qName, qLength, qStart, qEnd, qStrand, tName, tLength, tStart, tEnd, tStrand, score, numMatch, numMismatch, numIns, numDel, mapQV, qAlignedSeq, matchPattern, tAlignedSeq = [x for x in line.split(" ") if len(x)>0]
            tcov = (int(tEnd)-int(tStart))/float(tLength)
            if tcov < 0.9:
                continue
            ident = float(numMatch)/(int(qEnd)-int(qStart))
            spec, marker = tName.split("/")[0].rsplit("_", 1)
            assert marker == wildcards.ref
            mData = (spec, ident, float(score))
            try:
                data[qName].append(mData)
            except KeyError:
                data[qName] = [mData]
            
        with open(output[0], "w") as out:
            for qName, entryList in data.items():
                if len(entryList) == 1:
                    out.write(qName)
                    out.write("\t%s\t%s\t%s\n" % entryList[0])
                    continue
                entryList.sort(key=lambda x: x[1])
                best = entryList[0]
                for entry in entryList[1:]:
                    iDiff = best[1] - entry[1]
                    if iDiff > 0.01:
                        #this entry is 1% points less similar than the best: stop here
                        break
                    else:
                        out.write(qName)
                        out.write("\t%s\t%s\t%s\n" % entry)
                
                
rule compareCluCls:
    input: reads="clusters/{sample}_cluster.fasta", its="mapping/cluMatches/match_{sample}_ITS.tsv", ssu="mapping/cluMatches/match_{sample}_SSU.tsv", lsu="mapping/cluMatches/match_{sample}_LSU.tsv"
    output: "mapping/{sample}_cls_clusters.tsv"
    run:
        its = {}
        for line in open(input.its):
            read, cls, ident, scr = line.strip().split("\t")
            try:
                its[read].append(cls)
            except KeyError:
                its[read] = [cls]
        ssu = {}
        for line in open(input.ssu):
            read, cls, ident, scr = line.strip().split("\t")
            try:
                ssu[read].append(cls)
            except KeyError:
                ssu[read] = [cls]
        lsu = {}
        for line in open(input.lsu):
            read, cls, ident, scr = line.strip().split("\t")
            try:
                lsu[read].append(cls)
            except KeyError:
                lsu[read] = [cls]
        
        unclear=0
        unknown=0
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input.reads), "fasta"):
                comb = set(its.get(rec.id,[])) & set(ssu.get(rec.id,[])) & set(lsu.get(rec.id,[]))
                if len(comb) == 1:
                    cls = comb.pop()
                elif len(comb) == 0:
                    cls = "unknown"
                    unknown+=1
                else:
                    cls = "|".join(comb)
                    unclear+=1
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (samples[wildcards.sample],
                                                        rec.id, 
                                                        "|".join(set(ssu.get(rec.id,["unknown"]))), 
                                                        "|".join(set(its.get(rec.id,["unknown"]))), 
                                                        "|".join(set(lsu.get(rec.id,["unknown"]))),
                                                        cls))
        print("%s: %i unknown, %i unclear" % (wildcards.sample, unknown, unclear))
        

rule removeChimeraRef:
    input: seqs="primers/{sample}_primer.fasta", ref="fullRef/isolateSeqs.fasta"
    output: fasta="chimeraRef/{sample}.nochimera.fasta", tsv="chimeraRef/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_refChimera.log"
    shell:
        "%(vsearch)s --uchime_ref {input.seqs} --db {input.ref} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

def fullRefInput(wildcards):
    return ["consensus/%s_consensus.fasta" % s for s in isolates[wildcards.spec]]

rule getFullRef:
    input: fullRefInput
    output: "fullRef/{spec}.fasta"
    params: minSize=10
    log: "logs/fullRef_{spec}.log"
    run:
        seq=[]
        for r in range(len(input)):
            seqList = []
            for rec in SeqIO.parse(open(input[r]), "fasta"):
                size = int(rec.id.strip(";").rsplit(";", 1)[-1].split("=")[1])
                if size >= params.minSize:
                    seqList.append(rec)
            assert len(seqList) == 1
            seq.append(seqList[0])
        if len(input) > 1:
            with open(log[0], "w") as logFile:
                try:
                    assert seq[0].seq == seq[1].seq
                except AssertionError:
                    logFile.write("Different sequences from replicates:\n\n")
                    logFile.write(seq[0].format("fasta"))
                    logFile.write(seq[1].format("fasta"))
        seq[0].id = "isolate_%s" % wildcards.spec
        open(output[0], "w").write(seq[0].format("fasta"))

rule collectFullRef:
    input: expand("fullRef/{spec}.fasta", spec=isolates.keys())
    output: "fullRef/isolateSeqs.fasta"
    shell:
        "cat {input} > {output}"

rule fullRawMapping:
    input: reads="raw/{sample}.fastq", ref="fullRef/isolateSeqs.fasta"
    output: m5="mapping/fullMapping/{sample}_raw_vs_isolates.m5"
    threads: 6
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule fullMapping:
    input: reads="chimeraRef/{sample}.nochimera.fasta", ref="fullRef/isolateSeqs.fasta"
    output: m5="mapping/fullMapping/{sample}_filtered_vs_isolates.m5"
    threads: 6
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule getFullRawCls:
    input: m5="mapping/fullMapping/{sample}_raw_vs_isolates.m5", fastq="raw/{sample}.fastq"
    output: matches="mapping/fullMatches/match_{sample}_raw_isolates.tsv", assign="mapping/assignment/{sample}_raw_assignments.tsv"
    run:
        data={}

        for line in open(input.m5):
            qName, qLength, qStart, qEnd, qStrand, tName, tLength, tStart, tEnd, tStrand, score, numMatch, numMismatch, numIns, numDel, mapQV, qAlignedSeq, matchPattern, tAlignedSeq = [x for x in line.split(" ") if len(x)>0]
            qcov = (int(qEnd)-int(qStart))/float(qLength)
            if qcov < 0.9:
                continue
            ident = float(numMatch)/(int(qEnd)-int(qStart))
            spec = tName.split("_")[1]
            mData = (str(ident), score, str(int(qEnd)-int(qStart)), numMismatch, numIns, numDel)
            rId = qName.rsplit("/", 1)[0]
            try:
                qData = data[rId]
            except KeyError:
                data[rId] = {spec: mData}
            else:
                if not spec in qData or mData[0] > qData[spec][0]:
                    qData[spec] = mData
        with open(output.matches, "w") as out:
            for rId, entryDict in data.items():
                for spec, entry in entryDict.items():
                    out.write("%s\t%s\t%s\n" % (rId, spec, "\t".join(entry)))
        with open(output.assign, "w") as out:
            for rec in SeqIO.parse(open(input.fastq), "fastq"):
                rId = rec.id.rsplit("/", 1)[0]
                if rId in data:
                    matchDict = data[rId]
                    best = sorted(matchDict.items(), key=lambda x: int(x[1][1]))[0]
                    out.write("%s\t%s\t%s\n" % (rId, best[0], "\t".join(best[1][2:6])))
                else:
                    out.write("%s\tUNKNOWN\tNA\tNA\tNA\tNA\n" % rId)

rule getFullCls:
    input: m5="mapping/fullMapping/{sample}_filtered_vs_isolates.m5", chimera="chimeraRef/{sample}.chimeraReport.tsv"
    output: matches="mapping/fullMatches/match_{sample}_filtered_isolates.tsv", assign="mapping/assignment/{sample}_filtered_assignments.tsv"
    run:
        data={}

        for line in open(input.m5):
            qName, qLength, qStart, qEnd, qStrand, tName, tLength, tStart, tEnd, tStrand, score, numMatch, numMismatch, numIns, numDel, mapQV, qAlignedSeq, matchPattern, tAlignedSeq = [x for x in line.split(" ") if len(x)>0]
            qcov = (int(qEnd)-int(qStart))/float(qLength)
            if qcov < 0.9:
                continue
            ident = float(numMatch)/(int(qEnd)-int(qStart))
            spec = tName.split("_")[1]
            mData = (str(ident), score, str(int(qEnd)-int(qStart)), numMismatch, numIns, numDel)
            rId = qName.rsplit("/", 1)[0]
            try:
                qData = data[rId]
            except KeyError:
                data[rId] = {spec: mData}
            else:
                if not spec in qData or mData[0] > qData[spec][0]:
                    qData[spec] = mData
            
        with open(output.matches, "w") as out:
            for rId, entryDict in data.items():
                for spec, entry in entryDict.items():
                    out.write("%s\t%s\t%s\n" % (rId, spec, "\t".join(entry)))
        with open(output.assign, "w") as out:
            for line in open(input.chimera):
                arr = line.strip().split("\t")
                rId = arr[1]
                if arr[-1] != "N":
                    out.write("%s\tCHIMERA\tNA\tNA\tNA\tNA\n" % rId)
                elif rId in data:
                    matchDict = data[rId]
                    best = sorted(matchDict.items(), key=lambda x: int(x[1][1]))[0]
                    out.write("%s\t%s\t%s\n" % (rId, best[0], "\t".join(best[1][2:6])))
                else:
                    out.write("%s\tUNKNOWN\tNA\tNA\tNA\tNA\n" % rId)

rule collectFullCls:
    input: expand("mapping/assignment/{sample}_{{stage}}_assignments.tsv", sample=samples.keys())
    output: "mapping/all_{stage}_assignments.tsv"
    run:
        with open(output[0], "w") as out:
            for inputFile in input:
                sample = inputFile.rsplit("/", 1)[-1].split("_", 1)[0]
                for line in open(inputFile):
                    out.write("%s\t%s\t%s" % (sample, samples[sample], line))

rule plotErrorRate:
    input: "mapping/all_{stage}_assignments.tsv"
    output: "mapping/{stage}ErrorRates.svg"
    run:
        R("""
        library(ggplot2)
        library(reshape2)

        d=read.table("{input}")

        colnames(d) = c("sampleId", "sampleName", "readID", "cls", "alnLen", "sub", "ins", "del")
        d$pSub = d$sub/d$alnLen
        d$pIns = d$ins/d$alnLen
        d$pDel = d$del/d$alnLen

        e=subset(d, sampleName!="mock" & cls!="CHIMERA" & cls!="UNKNOWN")
        p=melt(e, id.vars=c("sampleId", "sampleName", "readID", "cls"), measure.vars=c("pSub", "pIns", "pDel"))
        ggplot(p) + geom_boxplot(aes(x=sampleName, y=value, fill=variable)) + labs(x="isolate sample", y="error rate", fill="error type") + scale_fill_brewer(labels=c("substitutions", "insertions", "deletions"), type="div") + theme_bw()
        ggsave("{output}", width=16, height=10)
        """)

rule plotAssignment:
    input: "mapping/all_{stage}_assignments.tsv"
    output: "mapping/{stage}Assignments.svg"
    run:
        R("""
        library(ggplot2)
        
        d=read.table("{input}")
        colnames(d) = c("sampleId", "sampleName", "readID", "cls", "alnLen", "sub", "ins", "del")
        d$cls=factor(d$cls, levels=c("UNKNOWN", "CHIMERA", "MR", "LS", "ME", "PB", "Psp", "CL", "SC", "DT", "CR", "TR", "PC", "Csp", "CHY1", "UM", "EV", "CA"))
        
        mycols = adjustcolor(rainbow(16), red.f = 0.7, blue.f = 0.7, green.f = 0.7)
        mycols = c("grey", "black", mycols)
        ggplot(d) + geom_bar(aes(x=sampleName, fill=cls)) +scale_fill_manual(values=mycols) + labs(x="sample", y="read number", fill="classification") + theme_bw()
        ggsave("{output}", width=16, height=10)
        """)
