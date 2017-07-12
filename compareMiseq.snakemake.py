import pickle
import math

from Bio import SeqIO

from snakemake.utils import min_version, R

shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow network storage

configfile: "config.json"

#miseq sample names are sorted, such that the first is the same replicate as used for the PacBio
sampleNames = {"STN_w_l" : {"PacBio": ("Lib3-0075",), "MiSeq": ("46_S121", "45_S120", "47_S122")},
               "STN_w_p" : {"PacBio": ("Lib3-0034",), "MiSeq": ("50_S124", "39_S114", "54_S128")},
               "STN_s_l" : {"PacBio": ("Lib7-0075",), "MiSeq": ("136_S35", "128_S61", "168_S74")},
               "STN_s_p" : {"PacBio": ("Lib7-0034",), "MiSeq": ("98_S53", "91_S26", "144_S29")}
              }
pb2ms = {}
ms2pb = {}
for sampleDict in sampleNames.values():
    pb2ms[sampleDict["PacBio"][0]] = sampleDict["MiSeq"]
    for m in sampleDict["MiSeq"]:
        ms2pb[m] = sampleDict["PacBio"][0]

rule all:
    input: expand("miseqComp/allOtus{ident}.size.tsv", ident=[95, 97, 99]), expand( "miseqComp/taxonomy/all_{ident}otu_ITS.class.tsv", ident=[95, 97, 99])

rule collectMsReads: 
    input: sampleInfo="../mycolinkMetabarcoding2/readInfo/sample_R1.tsv", repseq="../mycolinkMetabarcoding2/readInfo/all.repseq.tsv", fasta="../mycolinkMetabarcoding2/itsx/all.ITS2.fasta"
    output: fasta="miseqComp/its2/{sample}_msITS2.fasta", size="miseqComp/its2/{sample}_msITS2.info.tsv"
    run:
        sampleInf = {}
        for line in open(input.sampleInfo):
            read, sample = line.strip().split("\t")
            if sample == sampleNames[wildcards.sample]["MiSeq"][0]:
                sampleInf[read] = sample
        repSeqList = {}
        for line in open(input.repseq):
            seq, rep = line.strip().split("\t")
            try:
                repSeqList[rep].append(seq)
            except KeyError:
                repSeqList[rep] = [seq]
        with open(output.fasta, "w") as out, open(output.size, "w") as sizeOut:
            for rec in SeqIO.parse(open(input.fasta), "fasta"):
                name = rec.id.split("|")[0]
                seqId, sizeStr = name.strip(";").split(";")
                n=0
                for seq in repSeqList[seqId]:
                    if seq in sampleInf:
                        n+=1
                if n!=0:
                    rec.id = "%s;size=%i;" % (seqId, n)
                    out.write(rec.format("fasta"))
                    sizeOut.write("%s\tMiSeq\t%s\t%i\n" % (seqId, wildcards.sample, n))

def pbReadsInput(wildcards):
    return "itsx/%s.ITS2.fasta" % sampleNames[wildcards.sample]["PacBio"][0]

rule collectPbReads:
    input: pbReadsInput
    output: fasta="miseqComp/its2/{sample}_pbITS2.fasta", size="miseqComp/its2/{sample}_pbITS2.info.tsv"
    run:
        with open(output.fasta, "w") as out, open(output.size, "w") as sizeOut:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                seqId = rec.id.split("=", 1)[1].split(";", 1)[0]
                size = int(rec.id.split(";")[2].split("=")[1])
                rec.id = "%s_%s;size=%i;" % (wildcards.sample, seqId, size)
                out.write(rec.format("fasta"))
                sizeOut.write("%s_%s\tPacbio\t%s\t%i\n" % (wildcards.sample, seqId, wildcards.sample, size))

rule concatReads:
    input: expand("miseqComp/its2/{sample}_{tech}ITS2.fasta", sample=sampleNames.keys(), tech=["ms", "pb"])
    output: "miseqComp/its2/all.fasta"
    shell:
        "cat {input} > {output}"

rule concatInfo:
    input: expand("miseqComp/its2/{sample}_{tech}ITS2.info.tsv", sample=sampleNames.keys(), tech=["ms", "pb"])
    output: "miseqComp/its2/all.info.tsv"
    shell:
        "cat {input} > {output}"

rule cluster:
    input: "miseqComp/its2/all.fasta"
    output: fasta="miseqComp/otus/all.otus{ident}.fasta", uc="miseqComp/otus/all.otus{ident}.uc.txt"
    log: "miseqComp/logs/all_cluster{ident}.log"
    threads: 6
    shell:
        "%(vsearch)s --cluster_size {input} --relabel otu --sizein --sizeout --iddef 0 --id 0.{wildcards.ident} --minsl 0.9 --centroids {output.fasta} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

rule otuLength:
    input: "miseqComp/otus/all.otus{ident}.fasta"
    output: "miseqComp/otus/all.otus{ident}.length.tsv"
    run:
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                out.write("%s\t%i\n" % (rec.id.split(";")[0], len(rec)))

rule clusterInfo:
    input: uc="miseqComp/otus/all.otus{ident}.uc.txt", info="miseqComp/its2/all.info.tsv", cls="miseqComp/taxonomy/all_{ident}otu_ITS.class.tsv"
    output: "miseqComp/allOtus{ident}.size.tsv"
    run:
        info={}
        for line in open(input.info):
            seqId, machine, sample, size = line.strip().split("\t")
            info[seqId] = (machine, sample, size)
        clusters = {}
        for line in open(input.uc):
            if line[0] == "C":
                continue
            elif line[0] in "SH":
                arr = line.strip().split("\t")
                query, target = arr[-2:]
                otu = "otu%i" % (int(arr[1])+1)
                if target == "*":
                    #new cluster
                    clusters[otu] = [query]
                else:
                    clusters[otu].append(query)
            else:
                raise ValueError()
        cls={}
        for line in open(input.cls):
            otu, tCls = line.strip().split("\t")
            cls[otu] = tCls
        with open(output[0], "w") as out:
            for otu, seqList in clusters.items():
                otuSize = {s: {"MiSeq": 0, "Pacbio": 0} for s in sampleNames.keys()}
                for seq in seqList:
                    machine, sample, size = info[seq.split(";")[0]]
                    otuSize[sample][machine] += int(size)
                for tSample in sampleNames.keys():
                    ms = otuSize[tSample]["MiSeq"]
                    pb = otuSize[tSample]["Pacbio"]
                    phyl = "unknown"
                    try:
                        phyl=cls[otu].split(";")[1]
                    except IndexError:
                        pass
                    out.write("%s\t%s\t%s\t%s\t%i\n" % (otu, phyl, tSample, ms, pb))

rule plotOverlap:
    input: "miseqComp/allOtus{ident}.size.tsv"
    output: all="miseqComp/otu{ident}_overlap_all.svg", top="miseqComp/otu{ident}_overlap_top.svg"
    params: top=100
    run:
        R("""
        library(ggplot2)
        library(VennDiagram)

        d=read.table("{input}")
        colnames(d) = c("oId", "cls", "sample", "miseq", "pacbio")

        a=aggregate(cbind(miseq, pacbio)~oId+cls, d, sum)


        sum(a$miseq>0 & a$pacbio>0)
        # 530
        sum(a$miseq>0 & a$pacbio==0)
        # 1385
        sum(a$miseq==0 & a$pacbio>0)
        # 820

        svg("{output.all}")
        draw.pairwise.venn(sum(a$miseq>0), sum(a$pacbio>0), sum(a$miseq>0 & a$pacbio>0), c("Miseq", "PacBio"))
        dev.off()
        
        o=a[order(a$miseq+a$pacbio, decreasing=T),]
        tx=o[1:{params.top},]
        svg("{output.top}")
        draw.pairwise.venn(sum(tx$miseq>0), sum(tx$pacbio>0), sum(tx$miseq>0 & tx$pacbio>0), c("Miseq", "PacBio"))
        dev.off()
        
        """)

rule alignToUnite:
    input: clu="miseqComp/otus/all.otus{ident}.fasta", db="%(dbFolder)s/UNITE_%(uniteVersion)s.index.lambda" % config
    output: "miseqComp/lambda/all.{ident}otu_vs_UNITE.m8"
    log: "miseqComp/logs/all_{ident}otu_lambda.log"
    threads: 3
    resources: lamDb=1
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -p blastn -t {threads} &> {log}" % config

rule classifyITS:
    input: clu="miseqComp/otus/all.otus{ident}.fasta", lam="miseqComp/lambda/all.{ident}otu_vs_UNITE.m8", tax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
    output: "miseqComp/taxonomy/all_{ident}otu_ITS.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "miseqComp/logs/all_{ident}otuClass.log", "miseqComp/logs/all_{ident}otu_itsTax.log"
    run:
        taxDict = {}
        for line in open(input.tax):
            sh, tax = line.strip().split("\t")
            taxDict[sh] = tax
        logOut = open(log[0], "w")
        logTax = open(log[1], "w")
        classifi = {}
        seqLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
        for line in open(input.lam, encoding="latin-1"):
            total +=1
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.strip().split("\t")
            readId = qseqid.split(";")[0]
            if float(evalue) > params.maxE:
                evalueFilter += 1
                continue
            if float(pident) < params.minIdent:
                identFilter +=1
                continue
            mLen = min(int(qlen), int(slen))
            if float(length)/mLen*100 < params.minCov:
                covFilter += 1
                continue
            linStr = taxDict[sseqid]
            try:
                classifi[readId].append((linStr, float(bitscore)))
            except KeyError:
                classifi[readId] = [(linStr, float(bitscore))]
        logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
        logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, params.maxE))
        logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, params.minIdent))
        logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, params.minCov))
        topPerc = params.topPerc/100.0
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input.clu), "fasta"):
                otu = rec.id.split(";", 1)[0]
                if otu in classifi:
                    hits = classifi[otu]
                    sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                    cutoff = 0
                    while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                        cutoff += 1
                    goodHits = [hit[0] for hit in sortedHits[:cutoff]]
                    for h in goodHits:
                        logTax.write("%s\t%s\n" % (otu, h))
                    lineage = lca(goodHits, params.stringency)
                    out.write("%s\t%s\n" % (otu, lineage))
                else:
                    out.write("%s\tunknown\n" % (otu))
        try:
            logOut.close()
        except:
            pass
        try:
            logTax.close()
        except:
            pass

#rule generateMiSeqFile:
#    input: count="../mycolinkMetabarcoding2/swarm/{sampleId}.ITS2.otus.out", seqs="../mycolinkMetabarcoding2/swarm/all.ITS2.otus.fasta"
#    output: fasta="miseqOtus/{sampleId}.otus.fasta", size="miseqOtus/{sampleId}.otusize.pic", abund="miseqOtus/{sampleId}.otuAbund.pic"
#    run:
#        otuSeq = {}
#        otuSize = {}
#        for rec in SeqIO.parse(open(input.seqs), "fasta"):
#            otuId, sizeStr = rec.id.strip(";").split(";")
#            otuSeq[otuId] = rec
#            otuSize[otuId] = int(sizeStr.split("=")[1])
#        pickle.dump(otuSize, open(output.size, "wb"))
#        pickle.dump(otuAbund, open(output.abund, "wb"))
#        with open(output.fasta, "w") as out:
#            for line in open(input.count):
#                otuId, otuAbund = line.strip().split("\t")
#                seqRec = otuSeq[otuId]
#                seqRec.id = otuId
#                out.write(seqRec)

#rule createIndex:
#    input: "otus/{sampleId}_97otus.fasta"
#    output: touch("otus/{sampleId}_97otus.fasta.lambdaIndexCreated")
#    threads: 6
#    shell:
#        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config

#rule alignMiseqToPacBio:
#    input: miseq="../mycolinkMetabarcoding2/swarm/all.ITS2.otus.fasta", db="otus/{sampleId}_97otus.fasta", dbFlag="otus/{sampleId}_97otus.fasta.lambdaIndexCreated"
#    output: "miseqLambda/miseqOtu_vs_{sampleId}.m8"
#    log: "logs/miseqOtu_vs_{sampleId}_lambda.log"
#    threads: 3
#    shell:
#        "%(lambdaFolder)s/lambda -q {input.miseq} -d {input.db} -o {output} -p blastn -t {threads} &> {log}" % config

#rule assignMiseqToPacbio:
#    input: lam="miseqLambda/miseqOtu_vs_{pbSampleId}.m8", miseq="../mycolinkMetabarcoding2/swarm/all.ITS2.otus.fasta"
#    output: "miseqLambda/miseqOtu_vs_{pbSampleId}.assignment.tsv"
#    params: minCov=0.9, minIdent=97.0
#    run:
#        seqNr = 0
#        seqLength = {}
#        for rec in SeqIO.parse(open(input.miseq), "fasta"):
#            seqNr += 1
#            seqLength[rec.id] = len(rec)
#        matches = {}
#        covFilter = 0
#        for line in open(input[0]):
#            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
#            if float(pident) < params.minIdent:
#                continue
#            if float(length)/seqLength[qseqid] < params.minCov:
#                covFilter += 1
#                continue
#            try:
#                matches[qseqid].append((sseqid, pident))
#            except KeyError:
#                matches[qseqid] = [(sseqid, pident)]
#        with open(output[0], "w") as out:
#            for query, tMatches in matches.items():
#                sortedMatches = sorted(tMatches, key=lambda x: x[1], reverse=True)
#                out.write("%s\t%i\t%s\t%s\n" % (query, len(tMatches), 
#                                              ",".join([m[1] for m in sortedMatches]),
#                                              ",".join([m[0] for m in sortedMatches])))

#rule perRepAssign:
#    input: ass="miseqLambda/miseqOtu_vs_{pbSampleId}.assignment.tsv", msSize="../mycolinkMetabarcoding2/swarm/{msSampleId}.ITS2.otus.out" 
#    output: ass="miseqComp/{msSampleId}_vs_{pbSampleId}.assignment.tsv"
#    params: minSize=2
#    run:
#        msAbund = {}
#        for line in open(input.msSize):
#            otu, sizeStr = line.strip().split("\t")
#            size = int(sizeStr)
#            if size>=params.minSize:
#                msAbund[otu] = size
#        with open(output.ass, "w") as out:
#            for line in open(input.ass):
#                msOtu, mCount, mIdentList, pbOtusList = line.strip().split("\t")
#                msOtuId, msOtuSizeStr = msOtu.strip(";").split(";")
#                if msOtuId in msAbund:
#                    for pbOtuId in pbOtusList.split(","):
#                        out.write("%s\t%s\n" % (msOtuId, pbOtuId))
##                    best = pbOtusList.split(",")[0]
##                    out.write("%s\t%s\n" % (msOtuId, best))

#rule msAbund:
#    input: msSize="../mycolinkMetabarcoding2/swarm/{msSampleId}.ITS2.otus.out" 
#    output: abund="miseqComp/{msSampleId}.msOtuInfo.tsv"
#    run:
#        msAbund = {}
#        for line in open(input.msSize):
#            otu, sizeStr = line.strip().split("\t")
#            size = int(sizeStr)
#            if size>1:
#                msAbund[otu] = size
#        totalAbund = sum(msAbund.values())
#        with open(output.abund, "w") as out:
#            for otu, abund in msAbund.items():
#                out.write("%s:%s\t%f\tMiSeq\n" % (wildcards.msSampleId, otu, abund))

#rule pbAbund:
#    input: "otus/{pbSampleId}_97otus.size.tsv"
#    output: "miseqComp/{pbSampleId}.pbOtuInfo.tsv"
#    run:
#        pbAbund = {}
#        for line in open(input[0]):
#            otu, size = line.strip().split("\t")
#            pbAbund[otu] = int(size)
#        totalAbund = sum(pbAbund.values())
#        with open(output[0], "w") as out:
#            for otu, abund in pbAbund.items():
#                out.write("%s\t%f\tPacBio\n" % (otu, abund))

#def sampleAssInput(wildcards):
#    sample = sampleNames[wildcards.sample]
#    pbSampleId = sample["PacBio"][0]
#    msSampleIdList = sample["MiSeq"]
#    assFiles = ["miseqComp/%s_vs_%s.assignment.tsv" % (msSampleId, pbSampleId) for msSampleId in msSampleIdList]
#    return assFiles

#rule perSampleAssign:
#    input: sampleAssInput
#    output:"miseqComp/{sample}.assign.tsv"
#    run:
#        with open(output[0], "w") as out:
#            msNodes = {}
#            for inputFileName in input:
#                msSampleId, pbSampleId = inputFileName.split("/")[-1].split(".")[0].split("_vs_")
#                for line in open(inputFileName):
#                    msOtuId, pbOtuId = line.strip().split("\t")
#                    tNode = "%s:%s" % (msSampleId, msOtuId)
#                    out.write("%s\t%s/ITS\n" % (tNode, pbOtuId.split("|")[0]))
#                    #if we have never seen this miseq OTU ID before
#                    if msOtuId not in msNodes:
#                        msNodes[msOtuId] = [tNode]
#                    #if we have seen this miseq OTU ID before, but in a different replicate
#                    elif tNode not in msNodes[msOtuId]:
#                        for otherNode in msNodes[msOtuId]:
#                            out.write("%s\t%s\n" % (tNode, otherNode))
#                        msNodes[msOtuId].append(tNode)
#                        

#def sampleAbundInput(wildcards):
#    sample = sampleNames[wildcards.sample]
#    pbSampleId = sample["PacBio"][0]
#    msSampleIdList = sample["MiSeq"]
#    return ["miseqComp/%s.msOtuInfo.tsv" % msSampleId for msSampleId in msSampleIdList] + ["miseqComp/%s.pbOtuInfo.tsv" % pbSampleId]

#rule allAbund:
#    input: sampleAbundInput
#    output: "miseqComp/{sample}.otuInfo.tsv"
#    shell:
#        "cat {input} > {output}"

#rule sampleComp:
#    input: assign="miseqComp/{sample}.assign.tsv", info="miseqComp/{sample}.otuInfo.tsv"
#    output: ms="miseqComp/{sample}.msComp.tsv", pb="miseqComp/{sample}.pbComp.tsv"
#    run:
#        aInfo = {}
#        tInfo = {}
#        msSamples = []
#        for line in open(input.info):
#            otuId, abund, oType = line.strip().split("\t")
#            aInfo[otuId] = float(abund)
#            tInfo[otuId] = oType
#            if oType == "MiSeq":
#                msSample = otuId.split(":")[0]
#                if msSample not in msSamples:
#                    msSamples.append(msSample)
#        pbAssign = {}
#        msAssign = {}
#        for line in open(input.assign):
#            msOtu, pbOtu = line.strip().split("\t")
#            try:
#                pbAssign[pbOtu].append(msOtu)
#            except KeyError:
#                pbAssign[pbOtu] = [msOtu]
#            try:
#                assert msOtu not in msAssign
#            except AssertionError:
#                print(msOtu)
#                raise
#            msAssign[msOtu] = pbOtu
#        with open(output.ms, "w") as msOut, open(output.pb, "w") as pbOut:
#            for otu, abund in aInfo.items():
#                if tInfo[otu] == "PacBio":
#                    pbAbund = aInfo[otu]
#                    try:
#                        msOtuList = pbAssign[otu]
#                    except KeyError:
#                        msOtus = "NA"
#                        msAbund = {s: 0 for s in msSamples}
#                    else:
#                        msOtus = ",".join(msOtuList)
#                        msAbund = {s: 0 for s in msSamples}
#                        for msOtu in msOtuList:
#                            msSample, msId = msOtu.split(":")
#                            msAbund[msSample] += aInfo[msOtu]
#                    pbOut.write("%s\t%f\t%s\t" % (otu, pbAbund, msOtus))
#                    for s in msSamples:
#                        pbOut.write("%f\t" % msAbund[s])
#                    pbOut.write("\n")
#                else:
#                    msAbund = aInfo[otu]
#                    try:
#                        pbOtu = msAssign[otu]
#                    except KeyError:
#                        pbOtu = "NA"
#                        pbAbund = "NA"
#                    else:
#                        pbAbund = aInfo[pbOtu]
#                    msOut.write("%s\t%f\t%s\t%s\n" % (otu, msAbund, pbOtu, str(pbAbund)))


rule sampleGraph_createDot:
    input: edges="miseqComp/{sample}.assign.tsv", nodes="miseqComp/{sample}.otuInfo.tsv"
    output: "miseqComp/graph/{sample}.dot"
    params: maxNodeNr=20
    run:
        nodeId = {}
        usedNode = {}
        pbNodes = {}
        msNodes = {}
        pbNr=0
        msNr=0
        for line in open(input.nodes):
            name, abund, nType = line.strip().split("\t")
            if nType == "PacBio":
                pbNodes[name] = float(abund)
                nodeId[name] = "PB%i" % pbNr
                pbNr += 1
            else:
                msSampleId, msOtuId = name.split(":")
                try:
                    msNodes[msSampleId][name] = float(abund)
                except KeyError:
                    msNodes[msSampleId] = {name: float(abund)}
                nodeId[name] = "MS%i" % msNr
                msNr += 1
        dot='graph {{\nnode [style=filled, fontsize=4, fixedsize=true, shape="box", label=""];\n'\
            'edge [penwidth=2];\n'
        maxAbund = 0
        for rNr, msSampleId in enumerate(msNodes):
            maxAbund = max(maxAbund, max(msNodes[msSampleId].values()))
        scale = maxAbund/5
        for rNr, msSampleId in enumerate(msNodes):
            rank=0
            for otu, abund in sorted(msNodes[msSampleId].items(), key=lambda x: x[1], reverse=True):
                dot+='node [color=blue, height=%f, width=1, pos="%i,%i!"] %s;\n' % (abund/scale, rank*2, (rNr+1)*6, nodeId[otu])
                rank+=1
                usedNode[otu] = None
                if rank > params.maxNodeNr:
                    break
        rank=0
        for otu, abund in sorted(pbNodes.items(), key=lambda x: x[1], reverse=True):
            dot += 'node [color=red, height=%f, width=1, pos="%i,0!"] %s;\n' % (abund/scale, rank*2, nodeId[otu])
            rank += 1
            usedNode[otu] = None
            if rank > params.maxNodeNr:
                break
        for line in open(input.edges):
            start, end = line.strip().split("\t")
            if not start in usedNode or not end in usedNode:
                continue
            dot += '%s -- %s;\n' % (nodeId[start], nodeId[end])
        dot+="}}"
        open(output[0], "w").write(dot)
        
rule sampleGraph_graphiz:
    input: "miseqComp/graph/{sample}.dot"
    output: "miseqComp/graph/{sample}.svg"
    shell:
        "neato -Tsvg {input} -o {output}"



## HELPER FUNCTIONS
def lca(lineageStrings, stringency=1.0, 
        unidentified=["unidentified", "unclassified", "unknown"],
        ignoreIncertaeSedis=True, sizes=None):
    lineage = []
    mLineages = []
    #remove bootstrap values ("(100)", "(75)", etc.) if any
    for mLin in [l.strip(";").split(";") for l in lineageStrings]:
        mLineages.append([])
        for entry in mLin:
             mLineages[-1].append(entry.split("(")[0])
    maxLinLen = max([len(m) for m in mLineages])
    active = [True]*len(mLineages)
    for i in range(maxLinLen):
        total = 0.0
        counts = {}
        for m, memberLin in enumerate(mLineages):
            if not active[m]:
                continue #ignore lineages that were deactivated further up in the tree
            if len(memberLin) <= i:
                active[m] = False
                continue #ignore lineages that are not this long
            name = memberLin[i].split("__")[-1]
            if name in unidentified:
                continue # ignoring unidentified entrys
            if ignoreIncertaeSedis and name.startswith("Incertae"):
                continue # ignoring Incertae sedis entries.
                         # NOTE: this will mean lineages end at the first Incerta sedis
            if sizes is None:
                tSize = 1
            else:
                tSize = sizes[m]
            total += tSize
            try:
                counts[memberLin[i]] += tSize
            except KeyError:
                counts[memberLin[i]] = tSize
        if not counts:
            #no valid lineage entrys found in this level
            break
        most=sorted(counts.items(), key=lambda x: x[1], reverse=True)[0]
        different = total - most[1]
        #accept the lineage entry if its proportion of all (valid) classifications higher than stringency setting
        if different/total <= (1.0-stringency):
            lineage.append(most[0]) #add the most apearing entry to the new lineage
            #deactivate all lineages that were different at this level
            for m, memberLin in enumerate(mLineages):
                if active[m] and memberLin[i] != most[0]:
                    active[m] = False
        else:
            break
    if len(lineage) == 0:
        lineage = ["unknown"]
    return ";".join(lineage)
