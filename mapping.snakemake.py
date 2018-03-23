ref = ["LSU", "SSU", "ITS"]

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
    input: reads="%(inFolder)s/{sample}.fastq" % config, ref="fullRef/isolateSeqs.fasta"
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
    input: m5="mapping/fullMapping/{sample}_raw_vs_isolates.m5", fastq="%(inFolder)s/{sample}.fastq" % config
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

        env=c("STN_l_s", "STN_l_w", "STN_p_s", "STN_p_w")
        mock=c("emPCR", "c18i20", "c18i2", "c13i8", "c15i8", "c18i8", "c25i8", "c30i8")

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
        
        env=c("STN_l_s", "STN_l_w", "STN_p_s", "STN_p_w")
        
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
