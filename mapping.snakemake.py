ref = ["LSU", "SSU", "ITS"]

rule mapping:
    """deprecated"""
    input: reads="primers/{sample}_primer.fasta", ref="../PacBioMetabarcoding2/references/all_{ref}.fasta"
    output: m5="mapping/mapping/{sample}_vs_{ref}.m5"
    threads: 3
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule getCls:
    """deprecated"""
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
    """deprecated"""
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
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sampleName[wildcards.sample],
                                                        rec.id, 
                                                        "|".join(set(ssu.get(rec.id,["unknown"]))), 
                                                        "|".join(set(its.get(rec.id,["unknown"]))), 
                                                        "|".join(set(lsu.get(rec.id,["unknown"]))),
                                                        cls))
        print("%s: %i unknown, %i unclear" % (wildcards.sample, unknown, unclear))
        
rule concatCls:
    input: expand("mapping/{sample}_cls.tsv", sample=sampleName.keys())
    output: "mapping/combinedCls.tsv"
    shell:
        "cat {input} > {output}"
        
rule plot:
    """deprecated"""
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
    """deprecated"""
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
    """deprecated"""
    input: "mapping/byCls/{sample}-{cls}.fasta"
    output: "mapping/alnByCls/{sample}-{cls}.aln.fasta"
    log: "mapping/logs/byClsAln_{sample}-{cls}.log"
    threads: 4
    shell:
        "%(mafft)s --ep 1 --thread {threads} {input} > {output} 2> {log}" % config

        
rule clusterMapping:
    """deprecated"""
    input: reads="clusters/{sample}_cluster.fasta", ref="../PacBioMetabarcoding2/references/all_{ref}.fasta"
    output: m5="mapping/cluMapping/{sample}_vs_{ref}.m5"
    threads: 3
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule getCluCls:
    """deprecated"""
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
    """deprecated"""
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
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sampleName[wildcards.sample],
                                                        rec.id, 
                                                        "|".join(set(ssu.get(rec.id,["unknown"]))), 
                                                        "|".join(set(its.get(rec.id,["unknown"]))), 
                                                        "|".join(set(lsu.get(rec.id,["unknown"]))),
                                                        cls))
        print("%s: %i unknown, %i unclear" % (wildcards.sample, unknown, unclear))
        

#rule removeChimeraRef:
#    """Remove chimeras with vsearch reference based algorithm"""
#    input: seqs="primers/{sample}_primer.fasta", ref="fullRef/isolateSeqs.fasta"
#    output: fasta="chimeraRef/{sample}.nochimera.fasta", tsv="chimeraRef/{sample}.chimeraReport.tsv"
#    log: "logs/{sample}_refChimera.log"
#    shell:
#        "%(vsearch)s --uchime_ref {input.seqs} --db {input.ref} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

def fullRefInput(wildcards):
    """determine input data for fullRef rule according to isolate samples in the config"""
    return ["consensus/%s_consensus.fasta" % s for s in isolates[wildcards.spec]]

rule getFullRef:
    """Select one consensus per isolate sample"""
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
    """Collect all isolate sequence into one fasta file"""
    input: expand("fullRef/{spec}.fasta", spec=[s for s in isolates.keys() if s != "IF"])
    output: "fullRef/isolateSeqs.fasta"
    shell:
        "cat {input} > {output}"

rule fullRefITSx:
    """Use ITSx to find rRNA genes/regions in each isolate"""
    input: "fullRef/isolateSeqs.fasta"
    output: "fullRef/itsx/isolateSeqs.SSU.fasta", "fullRef/itsx/isolateSeqs.ITS1.fasta", "fullRef/itsx/isolateSeqs.5_8S.fasta", "fullRef/itsx/isolateSeqs.ITS2.fasta", "fullRef/itsx/isolateSeqs.LSU.fasta", "fullRef/itsx/isolateSeqs.summary.txt", "fullRef/itsx/isolateSeqs.positions.txt", "fullRef/itsx/isolateSeqs.full.fasta"
    threads: 6
    log: "fullRef/logs/isolateSeq_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o fullRef/itsx/isolateSeqs --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 500 -E 1e-4 2> {log}" % config

rule createFullRefAnnoation:
    """Create gff annotation file with found rRNA genes/regions"""
    input: pos="fullRef/itsx/isolateSeqs.positions.txt"
    output: anno="fullRef/isolateSeqs.gff3"
    run:
        typeDict = {"LSU": "large_subunit_rRNA",
                    "SSU": "small_subunit_rRNA",
                    "5.8S": "rRNA_5_8S"}
                    
        with open(output.anno, "w") as out:
            for l, line in enumerate(open(input.pos)):
                seqid, lenStr, regions = line.strip().split("\t", 2)
                for regStr in regions.split("\t"):
                    nameStr, rangeStr = regStr.split(": ")
                    start, end = rangeStr.split("-")
                    name = nameStr.strip(":")
                    
                    if name in typeDict:
                        #seqid, source, type, start, end, score, strand, phase, attr
                        out.write("%s\tITSx\t%s\t%s\t%s\t.\t+\t.\t%s\n" % (seqid, typeDict[name], start, end, "ID=%i;Name=%s" % (l, name)))

rule fullRawMapping:
    """Map raw reads against isolate sequences with blasr"""
    input: reads="raw/{sample}.fastq", ref="fullRef/isolateSeqs.fasta"
    output: m5="mapping/fullMapping/{sample}_raw_vs_isolates.m5"
    threads: 6
    shell: 
        "/home/heeger/bin/blasr/blasr -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}"

rule fullMapping:
    """Map filtered reads against isolate sequences with blasr"""
    input: reads="refChimera/{sample}.nochimera.fasta", ref="fullRef/isolateSeqs.fasta"
    output: m5="mapping/fullMapping/{sample}_filtered_vs_isolates.m5"
    threads: 6
    shell: 
        "%(blasr)s -m 5 --bestn 50 --nproc {threads} --minPctSimilarity 90 --out {output.m5} {input.reads} {input.ref}" % config

rule getFullRawCls:
    """Create table with assignment (including chimera) for each raw read"""
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
    """Create table with assignment (including chimera) for each filtered read"""
    input: m5="mapping/fullMapping/{sample}_filtered_vs_isolates.m5", chimera="refChimera/{sample}.chimeraReport.tsv"
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
    """Collect mapping data from all samples into one table"""
    input: expand("mapping/assignment/{sample}_{{stage}}_assignments.tsv", sample=sampleName.keys())
    output: "mapping/all_{stage}_assignments.tsv"
    run:
        with open(output[0], "w") as out:
            for inputFile in input:
                sample = inputFile.rsplit("/", 1)[-1].split("_", 1)[0]
                for line in open(inputFile):
                    out.write("%s\t%s\t%s" % (sample, sampleName[sample], line))

rule plotErrorRate:
    """Plot error rate per type and per sample as a boxplot"""
    input: "mapping/all_{stage}_assignments.tsv"
    output: "mapping/{stage}ErrorRates.svg"
    run:
        R("""
        library(ggplot2)
        library(reshape2)

        env=c("DGW_l_s", "DGW_l_w", "DGW_p_s", "DGW_p_w", "FUKUSW_l_s", "FUKUSW_l_w", "FUKUSW_p_s", "GRB_l_s", "GRB_l_w", "GRB_p_s", "GRB_p_w", "STN_l_s", "STN_l_w", "STN_p_s", "STN_p_w")
        mock=paste("mock", c("emPCR", "i20c18", "i2c18", "i8c13", "i8c15", "i8c18", "i8c25", "i8c30"), sep="_")

        d=read.table("{input}")

        colnames(d) = c("sampleId", "sampleName", "readID", "cls", "alnLen", "sub", "ins", "del")
        d$pSub = d$sub/d$alnLen
        d$pIns = d$ins/d$alnLen
        d$pDel = d$del/d$alnLen

        e=subset(d, !(sampleName %in% mock) & !(sampleName %in% env) & cls!="CHIMERA" & cls!="UNKNOWN")
        p=melt(e, id.vars=c("sampleId", "sampleName", "readID", "cls"), measure.vars=c("pSub", "pIns", "pDel"))
        ggplot(p) + geom_boxplot(aes(x=sampleName, y=value, fill=variable)) + labs(x="isolate sample", y="error rate", fill="error type") + scale_fill_brewer(labels=c("substitutions", "insertions", "deletions"), type="div") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output}", width=16, height=10)
        """)

rule plotAssignment:
    """Plot occurence of each isoalte per sample as a heat map and composition of 
    each sample as a bar plot (for raw or filterd reads)"""
    input: "mapping/all_{stage}_assignments.tsv"
    output: ass="mapping/{stage}Assignments.svg", occ="mapping/{stage}Occurence.svg", mock="mapping/{stage}MockComp.svg"
    run:
        R("""
        library(ggplot2)
        
        env=c("DGW_l_s", "DGW_l_w", "DGW_p_s", "DGW_p_w", "FUKUSW_l_s", "FUKUSW_l_w", "FUKUSW_p_s", "GRB_l_s", "GRB_l_w", "GRB_p_s", "GRB_p_w", "STN_l_s", "STN_l_w", "STN_p_s", "STN_p_w")
        
        d=read.table("{input}")
        colnames(d) = c("sampleId", "sampleName", "readID", "cls", "alnLen", "sub", "ins", "del")
        d$cls=factor(d$cls, levels=c("UNKNOWN", "CHIMERA", "MR", "LS", "ME", "PB", "Psp", "CL", "SC", "DT", "CR", "TR", "PC", "Csp", "CHY1", "UM", "EV", "CA"))
        
        p=subset(d, !(sampleName %in% env))
        
        mycols = adjustcolor(rainbow(16), red.f = 0.7, blue.f = 0.7, green.f = 0.7)
        mycols = c("grey", "black", mycols)
        ggplot(p) + geom_bar(aes(x=sampleName, fill=cls)) +scale_fill_manual(values=mycols) + labs(x="sample", y="read number", fill="classification") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output.ass}", width=16, height=10)
        
        a=aggregate(readID~cls+sampleName, p, length)
        ggplot(a) + geom_tile(aes(sampleName, cls, fill=log10(readID))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x="sample", y="classification", fill="log10(read number)") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output.occ}", width=16, height=10)
        
        ggplot(p) + geom_bar(aes(x=sampleName, fill=cls), position="fill") +scale_fill_manual(values=mycols, drop=F) + labs(x="sample", y="read number", fill="classification") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output.mock}", width=16, height=10)
        """)
