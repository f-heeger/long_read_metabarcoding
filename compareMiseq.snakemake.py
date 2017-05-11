import pickle
import math

from Bio import SeqIO

from snakemake.utils import min_version, R

shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow network storage

configfile: "config.json"

sampleNames = {"STN_w_l" : {"PacBio": ("Lib3-0075",), "MiSeq": ("45_S120", "46_S121", "47_S122")},
               "STN_w_p" : {"PacBio": ("Lib3-0034",), "MiSeq": ("39_S114", "50_S124", "54_S128")},
               "STN_s_l" : {"PacBio": ("Lib7-0075",), "MiSeq": ("128_S61", "136_S35", "168_S74")},
               "STN_s_p" : {"PacBio": ("Lib7-0034",), "MiSeq": ("91_S26", "98_S53", "144_S29")}
              }
pb2ms = {}
for sampleDict in sampleNames.values():
    pb2ms[sampleDict["PacBio"][0]] = sampleDict["MiSeq"]

rule all:
    input: expand("miseqLambda/miseqOtu_vs_{sampleId}.assignment.tsv", sampleId=[s["PacBio"][0] for s in sampleNames.values()])

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

rule createIndex:
    input: "otus/{sampleId}_97otus.fasta"
    output: touch("otus/{sampleId}_97otus.fasta.lambdaIndexCreated")
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config

rule alignMiseqToPacBio:
    input: miseq="../mycolinkMetabarcoding2/swarm/all.ITS2.otus.fasta", db="otus/{sampleId}_97otus.fasta", dbFlag="otus/{sampleId}_97otus.fasta.lambdaIndexCreated"
    output: "miseqLambda/miseqOtu_vs_{sampleId}.m8"
    log: "logs/miseqOtu_vs_{sampleId}_lambda.log"
    threads: 3
    shell:
        "%(lambdaFolder)s/lambda -q {input.miseq} -d {input.db} -o {output} -p blastn -t {threads} &> {log}" % config

rule assignMiseqToPacbio:
    input: lam="miseqLambda/miseqOtu_vs_{pbSampleId}.m8", miseq="../mycolinkMetabarcoding2/swarm/all.ITS2.otus.fasta"
    output: "miseqLambda/miseqOtu_vs_{pbSampleId}.assignment.tsv"
    params: minCov=0.9, minIdent=97.0
    run:
        seqNr = 0
        seqLength = {}
        for rec in SeqIO.parse(open(input.miseq), "fasta"):
            seqNr += 1
            seqLength[rec.id] = len(rec)
        matches = {}
        covFilter = 0
        for line in open(input[0]):
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
            if float(pident) < params.minIdent:
                continue
            if float(length)/seqLength[qseqid] < params.minCov:
                covFilter += 1
                continue
            try:
                matches[qseqid].append((sseqid, pident))
            except KeyError:
                matches[qseqid] = [(sseqid, pident)]
        with open(output[0], "w") as out:
            for query, tMatches in matches.items():
                sortedMatches = sorted(tMatches, key=lambda x: x[1], reverse=True)
                out.write("%s\t%i\t%s\t%s\n" % (query, len(tMatches), 
                                              ",".join([m[1] for m in sortedMatches]),
                                              ",".join([m[0] for m in sortedMatches])))

rule perRepAssign:
    input: ass="miseqLambda/miseqOtu_vs_{pbSampleId}.assignment.tsv", msSize="../mycolinkMetabarcoding2/swarm/{msSampleId}.ITS2.otus.out" 
    output: ass="miseqComp/{msSampleId}_vs_{pbSampleId}.assignment.tsv"
    params: minSize=2
    run:
        msAbund = {}
        for line in open(input.msSize):
            otu, sizeStr = line.strip().split("\t")
            size = int(sizeStr)
            if size>=params.minSize:
                msAbund[otu] = size
        with open(output.ass, "w") as out:
            for line in open(input.ass):
                msOtu, mCount, mIdentList, pbOtusList = line.strip().split("\t")
                msOtuId, msOtuSizeStr = msOtu.strip(";").split(";")
                if msOtuId in msAbund:
                    for pbOtuId in pbOtusList.split(","):
                        out.write("%s\t%s\n" % (msOtuId, pbOtuId))
#                    best = pbOtusList.split(",")[0]
#                    out.write("%s\t%s\n" % (msOtuId, best))

rule msAbund:
    input: msSize="../mycolinkMetabarcoding2/swarm/{msSampleId}.ITS2.otus.out" 
    output: abund="miseqComp/{msSampleId}.msOtuInfo.tsv"
    run:
        msAbund = {}
        for line in open(input.msSize):
            otu, sizeStr = line.strip().split("\t")
            size = int(sizeStr)
            if size>1:
                msAbund[otu] = size
        totalAbund = sum(msAbund.values())
        with open(output.abund, "w") as out:
            for otu, abund in msAbund.items():
                out.write("%s:%s\t%f\tMiSeq\n" % (wildcards.msSampleId, otu, abund))

rule pbAbund:
    input: "otus/{pbSampleId}_97otus.size.tsv"
    output: "miseqComp/{pbSampleId}.pbOtuInfo.tsv"
    run:
        pbAbund = {}
        for line in open(input[0]):
            otu, size = line.strip().split("\t")
            pbAbund[otu] = int(size)
        totalAbund = sum(pbAbund.values())
        with open(output[0], "w") as out:
            for otu, abund in pbAbund.items():
                out.write("%s\t%f\tPacBio\n" % (otu, abund))

def sampleAssInput(wildcards):
    sample = sampleNames[wildcards.sample]
    pbSampleId = sample["PacBio"][0]
    msSampleIdList = sample["MiSeq"]
    assFiles = ["miseqComp/%s_vs_%s.assignment.tsv" % (msSampleId, pbSampleId) for msSampleId in msSampleIdList]
    return assFiles

rule perSampleAssign:
    input: sampleAssInput
    output:"miseqComp/{sample}.assign.tsv"
    run:
        with open(output[0], "w") as out:
            msNodes = {}
            for inputFileName in input:
                msSampleId, pbSampleId = inputFileName.split("/")[-1].split(".")[0].split("_vs_")
                for line in open(inputFileName):
                    msOtuId, pbOtuId = line.strip().split("\t")
                    tNode = "%s:%s" % (msSampleId, msOtuId)
                    out.write("%s\t%s/ITS\n" % (tNode, pbOtuId.split("|")[0]))
                    #if we have never seen this miseq OTU ID before
                    if msOtuId not in msNodes:
                        msNodes[msOtuId] = [tNode]
                    #if we have seen this miseq OTU ID before, but in a different replicate
                    elif tNode not in msNodes[msOtuId]:
                        for otherNode in msNodes[msOtuId]:
                            out.write("%s\t%s\n" % (tNode, otherNode))
                        msNodes[msOtuId].append(tNode)
                        

def sampleAbundInput(wildcards):
    sample = sampleNames[wildcards.sample]
    pbSampleId = sample["PacBio"][0]
    msSampleIdList = sample["MiSeq"]
    return ["miseqComp/%s.msOtuInfo.tsv" % msSampleId for msSampleId in msSampleIdList] + ["miseqComp/%s.pbOtuInfo.tsv" % pbSampleId]

rule allAbund:
    input: sampleAbundInput
    output: "miseqComp/{sample}.otuInfo.tsv"
    shell:
        "cat {input} > {output}"

rule sampleComp:
    input: assign="miseqComp/{sample}.assign.tsv", info="miseqComp/{sample}.otuInfo.tsv"
    output: ms="miseqComp/{sample}.msComp.tsv", pb="miseqComp/{sample}.pbComp.tsv"
    run:
        aInfo = {}
        tInfo = {}
        msSamples = []
        for line in open(input.info):
            otuId, abund, oType = line.strip().split("\t")
            aInfo[otuId] = float(abund)
            tInfo[otuId] = oType
            if oType == "MiSeq":
                msSample = otuId.split(":")[0]
                if msSample not in msSamples:
                    msSamples.append(msSample)
        pbAssign = {}
        msAssign = {}
        for line in open(input.assign):
            msOtu, pbOtu = line.strip().split("\t")
            try:
                pbAssign[pbOtu].append(msOtu)
            except KeyError:
                pbAssign[pbOtu] = [msOtu]
            try:
                assert msOtu not in msAssign
            except AssertionError:
                print(msOtu)
                raise
            msAssign[msOtu] = pbOtu
        with open(output.ms, "w") as msOut, open(output.pb, "w") as pbOut:
            for otu, abund in aInfo.items():
                if tInfo[otu] == "PacBio":
                    pbAbund = aInfo[otu]
                    try:
                        msOtuList = pbAssign[otu]
                    except KeyError:
                        msOtus = "NA"
                        msAbund = {s: 0 for s in msSamples}
                    else:
                        msOtus = ",".join(msOtuList)
                        msAbund = {s: 0 for s in msSamples}
                        for msOtu in msOtuList:
                            msSample, msId = msOtu.split(":")
                            msAbund[msSample] += aInfo[msOtu]
                    pbOut.write("%s\t%f\t%s\t" % (otu, pbAbund, msOtus))
                    for s in msSamples:
                        pbOut.write("%f\t" % msAbund[s])
                    pbOut.write("\n")
                else:
                    msAbund = aInfo[otu]
                    try:
                        pbOtu = msAssign[otu]
                    except KeyError:
                        pbOtu = "NA"
                        pbAbund = "NA"
                    else:
                        pbAbund = aInfo[pbOtu]
                    msOut.write("%s\t%f\t%s\t%s\n" % (otu, msAbund, pbOtu, str(pbAbund)))


rule sampleGraph_createDot:
    input: edges="miseqComp/{sample}.assign.tsv", nodes="miseqComp/{sample}.otuInfo.tsv"
    output: "miseqComp/{sample}.dot"
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
    input: "miseqComp/{sample}.dot"
    output: "miseqComp/{sample}.svg"
    shell:
        "neato -Tsvg {input} -o {output}"
