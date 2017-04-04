from snakemake.utils import min_version, R
from Bio import SeqIO

min_version("3.5.4")

shell.prefix("sleep 10; ")

configfile: "config.json"

samples = {"Lib4_0018": "mock",
           "Lib1_0009": "CA1",
           "Lib1_0027": "PC1",
           "Lib1_0056": "Csp1",
           "Lib1_0095": "CL1",
           "Lib2_0009": "UM1",
           "Lib2_0027": "PC2",
           "Lib2_0056": "ME1",
           "Lib3_0009": "SC1",
           "Lib3_0027": "LS1",
           "Lib3_0056": "Csp2",
           "Lib3_0095": "DT1",
           "Lib4_0009": "MR1",
           "Lib4_0027": "EV1",
           "Lib4_0056": "ME2",
           "Lib5_0009": "CA2",
           "Lib5_0027": "EV2",
           "Lib5_0095": "DT2",
           "Lib6_0009": "UM2",
           "Lib6_0027": "IF1",
           "Lib6_0095": "CL2",
           "Lib7_0009": "SC2",
           "Lib7_0027": "LS2",
           "Lib7_0056": "PB1",
           "Lib7_0095": "Psp1",
           "Lib8_0009": "MR2",
           "Lib8_0027": "IF2",
           "Lib8_0056": "PB2",
           "Lib8_0095": "Psp2",
           }

ref = ["LSU", "SSU", "ITS"]

rule all:
    input: "mapping/readPerSampleBarplot.pdf"

rule mapping:
    input: reads="primers/{sample}_minLen.fasta", ref="../PacBioMetabarcoding2/references/all_{ref}.fasta"
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
                
                
rule compareCls:
    input: reads="primers/{sample}_minLen.fasta", its="mapping/matches/match_{sample}_ITS.tsv", ssu="mapping/matches/match_{sample}_SSU.tsv", lsu="mapping/matches/match_{sample}_LSU.tsv"
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
    input: reads="primers/{sample}_minLen.fasta", cls="mapping/combinedCls.tsv"
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
    input: seqs="primers/{sample}_minLen.fasta", ref="isolateSeqs.fasta"
    output: fasta="chimeraRef/{sample}.nochimera.fasta", tsv="chimeraRef/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_refChimera.log"
    shell:
        "%(vsearch)s --uchime_ref {input.seqs} --db {input.ref} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule fullMapping:
    input: reads="chimeraRef/{sample}.nochimera.fasta", ref="isolateSeqs.fasta"
    output: m5="mapping/fullMapping/{sample}_vs_isolates.m5"
    threads: 6
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule getFullCls:
    input: m5="mapping/fullMapping/{sample}_vs_isolates.m5", chimera="chimeraRef/{sample}.chimeraReport.tsv"
    output: matches="mapping/fullMatches/match_{sample}_isolates.tsv", assign="mapping/assignment/{sample}_assignments.tsv"
    run:
        data={}

        for line in open(input.m5):
            qName, qLength, qStart, qEnd, qStrand, tName, tLength, tStart, tEnd, tStrand, score, numMatch, numMismatch, numIns, numDel, mapQV, qAlignedSeq, matchPattern, tAlignedSeq = [x for x in line.split(" ") if len(x)>0]
            qcov = (int(qEnd)-int(qStart))/float(qLength)
            if qcov < 0.9:
                continue
            ident = float(numMatch)/(int(qEnd)-int(qStart))
            spec = tName.split("_")[2]
            mData = (ident, float(score))
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
                    out.write("%s\t%s\t%s\t%s\t%i\n" % (rId, spec, entry[0], entry[1], len(entryDict)))
        with open(output.assign, "w") as out:
            for line in open(input.chimera):
                arr = line.strip().split("\t")
                rId = arr[1]
                if arr[-1] != "N":
                    out.write("%s\tCHIMERA\n" % rId)
                elif rId in data:
                    matchDict = data[rId]
                    best = sorted(matchDict.items(), key=lambda x: x[1][1])[0]
                    out.write("%s\t%s\n" % (rId, best[0]))
                else:
                    out.write("%s\tUNKNOWN\n" % rId)
