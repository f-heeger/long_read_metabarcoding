import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from snakemake.utils import min_version, R

import networkx as nx

shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow samba

#include: "mapping.snakemake.py"

configfile: "config.json"

samples = expand("Lib{nr}-0075", nr=range(1,9)) + expand("Lib{nr}-0034", nr=[1,2,3,5,6,7,8])
stechlin = ["Lib3-0075", "Lib3-0034", "Lib7-0075", "Lib7-0034"]
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
            "IF": ["Lib6-0027", "Lib8-0027"],
            "MR": ["Lib4-0009", "Lib8-0009"]
            }

#samples = []
#for lib, bcList in config["samples"].items():
#    for bc in bcList: 
#        samples.append("%s-%s" % (lib, bc))
#samples=["Lib4-0018"]
#print(sorted(samples))


####################################################################
# includes

include: "createDBs.snakemake.py"
include: "readProcessing.snakemake.py"

####################################################################

rule all:
#    input: "readNumbers.pdf", expand("chimera/{sample}.nochimera.fasta", sample=samples)   
    input: "taxonomy/all_97_comb.class.tsv", "all_clsComp_depth.svg", "all_clsComp_depth_fungi.svg", "all_clsComp_basic.svg", "taxonomy/Lib4-0018_97_combToCorr.class.tsv", "taxonomy/isolates_97_comb.class.tsv"
#    input: expand("primers/{sample}_primer.fasta", sample=samples)

rule concatItsxResult:
    input: expand("itsx/{sample}.{{marker}}.fasta", sample=samples)
    output: "catItsx/all.{marker}.fasta"
    shell:
        "cat {input} > {output}"

rule concatStechlin:
    input: expand("itsx/{sample}.{{marker}}.fasta", sample=stechlin)
    output: "catItsx/stechlin.{marker}.fasta"
    shell:
        "cat {input} > {output}"

def catIsolatesInput(wildcards):
    return ["itsx/%s.%s.fasta" % (s, wildcards.marker) for s in isolates[wildcards.spec]]

rule concatIsolates:
        input: catIsolatesInput
        output:"isolates/{spec}.{marker}.fasta"
        shell:
            "cat {input} > {output}"

def otuInput(wildcards):
    if wildcards.sample == "all":
        return "catItsx/all.full.fasta"
    elif wildcards.sample == "stechlin":
        return "catItsx/stechlin.full.fasta"
    elif wildcards.sample in isolates:
        return "isolates/%s.full.fasta" % wildcards.sample
    else:
        return "itsx/%s.full.fasta" % wildcards.sample

rule otuCluster:
    input: otuInput
    output: fasta="otus/{sample}_{ident}otus.fasta", uc="otus/{sample}_{ident}otus.uc.tsv"
    log: "logs/{sample}_{ident}otuClustering.log"
    threads: 3
    shell:
        "%(vsearch)s --cluster_size {input} --relabel otu --sizein --sizeout --iddef 0 --id 0.{wildcards.ident} --minsl 0.9 --centroids {output.fasta} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

def otuReadsInput(wildcards):
    if wildcards.sample == "all":
        return ["otus/all_%sotus.uc.tsv" % wildcards.ident] + expand("preclusters/{sample}_cluInfo.tsv", sample=samples)
    elif wildcards.sample == "stechlin":
        return ["otus/stechlin_%sotus.uc.tsv" % wildcards.ident] + expand("preclusters/{sample}_cluInfo.tsv", sample=stechlin)
    elif wildcards.sample in isolates:
        return ["otus/%s_%sotus.uc.tsv" % (wildcards.sample, wildcards.ident)] + expand("preclusters/{sample}_cluInfo.tsv", sample=isolates[wildcards.sample])
    else:
        return ["otus/%s_%sotus.uc.tsv" % (wildcards.sample, wildcards.ident), "preclusters/all_cluInfo.tsv"]

rule otuReads:
    input: otuReadsInput
    output: size="otus/{sample}_{ident}otus.size.tsv", info="otus/{sample}_{ident}otuReadInfo.pic", preInfo="otus/{sample}_{ident}otu_preClusterInfo.pic", repInfo="otus/{sample}_{ident}otu_repSeq.pic"
    run:
        preClusterReads={}
        for preInfo in input[1:]:
            for line in open(preInfo):
                clu, seq = line.strip().split("\t")
                try:
                    preClusterReads[clu].append(seq)
                except KeyError:
                    preClusterReads[clu] = [seq]
        clusterReads = {}
        otuPreClusters = {}
        repSeq = {}
        for line in open(input[0]):
            if line[0] == "C":
                arr = line.strip().split("\t")
                otu = "otu%i" % (int(arr[1])+1)
                seq = arr[-2]
                seqId = seq.split("=", 1)[1].split(";", 1)[0]
                repSeq[otu] = seqId
            elif line[0] in "SH":
                arr = line.strip().split("\t")
                seq = arr[-2]
                seqId = seq.split("=", 1)[1].split(";", 1)[0]
                otu = "otu%i" % (int(arr[1])+1)
                try:
                    otuPreClusters[otu].append(seq)
                except KeyError:
                    otuPreClusters[otu] = [seq]
                for read in preClusterReads[seqId]:
                    try:
                        clusterReads[otu].append(read)
                    except KeyError:
                        clusterReads[otu] = [read]
            else:
                raise ValueError("Unknown record type: %s" % arr[0])
        otuInf = {}
        with open(output.size, "w") as sOut:
            for otu, readList in clusterReads.items():
                for read in readList:
                    otuInf[read] = otu
                sOut.write("%s\t%i\n" % (otu, len(readList)))
        with open(output.info, "wb") as pOut:
            pickle.dump(otuInf, pOut)
        with open(output.preInfo, "wb") as eOut:
            pickle.dump(otuPreClusters, eOut)
        with open(output.repInfo, "wb") as rOut:
            pickle.dump(repSeq, rOut)


def transferOtusInput(wildcards):
    if wildcards.sample == "all":
        return ["otus/all_%sotu_repSeq.pic" % wildcards.ident, "catItsx/all.%s.fasta" % wildcards.marker]
    elif wildcards.sample == "stechlin":
        return ["otus/stechlin_%sotu_repSeq.pic" % wildcards.ident, "catItsx/stechlin.%s.fasta" % wildcards.marker]
    elif wildcards.sample in isolates:
        return ["otus/%s_%sotu_repSeq.pic" % (wildcards.sample, wildcards.ident), "isolates/%s.%s.fasta" % (wildcards.sample, wildcards.marker)]
    else:
        return ["otus/%s_%sotu_repSeq.pic" % (wildcards.sample, wildcards.ident), "itsx/%s.%s.fasta" % (wildcards.sample, wildcards.marker)]


rule transferOtus:
    input: transferOtusInput
    output: "otus/{sample}_{ident}otus_{marker}.fasta"
    run:
        repSeq = pickle.load(open(input[0], "rb"))
        rep2otu = dict(zip(repSeq.values(), repSeq.keys()))
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input[1]), "fasta"):
                readId = rec.id.split("=", 1)[1].split(";",1)[0]
                if readId in rep2otu:
                    rec.id = "%s/%s" % (rep2otu[readId], wildcards.marker)
                    out.write(rec.format("fasta"))

rule alignToUnite:
    input: clu="otus/{sample}_{ident}otus.fasta", db="%(dbFolder)s/UNITE_%(uniteVersion)s.index.lambda" % config
    output: "lambda/{sample}.{ident}otu_vs_UNITE.m8"
    log: "logs/{sample}_{ident}otu_lambda.log"
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -p blastn -t {threads} &> {log}" % config

rule classifyITS:
    input: lam="lambda/{sample}.{ident}otu_vs_UNITE.m8", tax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
    output: "taxonomy/{sample}_{ident}otu_ITS.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/{sample}_{ident}otuClass.log", "logs/{sample}_{ident}otu_itsTax.log"
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
            readId = qseqid.split("|")[0]
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
            for key, hits in classifi.items():
                sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                cutoff = 0
                while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                    cutoff += 1
                goodHits = [hit[0] for hit in sortedHits[:cutoff]]
                for h in goodHits:
                    logTax.write("%s\t%s\n" % (key, h))
                lineage = lca(goodHits, params.stringency)
                out.write("%s\t%s\n" % (key, lineage))
        try:
            logOut.close()
        except:
            pass
        try:
            logTax.close()
        except:
            pass

rule alignToSilva:
    input: clu="otus/{sample}_{ident}otus_SSU.fasta", db="%(dbFolder)s/silva_SSU_index.lambda" % config
    output: "lambda/{sample}.{ident}otu_SSU_vs_SILVA.m8"
    log: "logs/{sample}_{ident}otu_SSU_lambda.log"
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -nm 20000 -p blastn -t {threads} -b -2 -x 30 -as F &> {log}" % config

rule classifySSU:
    input: lam="lambda/{sample}.{ident}otu_SSU_vs_SILVA.m8", tax="%(dbFolder)s/SILVA_%(silvaVersion)s_SSU_tax.tsv" % config
    output: "taxonomy/{sample}_{ident}otu_SSU.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/{sample}_SSU_{ident}otuClass.log", "logs/{sample}_{ident}otu_SSU_tax.log", "logs/{sample}_{ident}otu_SSU_tiling.log"
    run:
        taxDict={}
        for line in open(input.tax):
            rId, tax = line.strip().split("\t")
            taxDict[rId] = tax
        classifi = {}
        seqLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
        hsp = {}
        sLen = {}
        qLen = {}
        for line in open(input.lam, encoding="latin-1"):
            total +=1
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.strip().split("\t")
            otuId = qseqid.split("/")[0]
            if float(evalue) > params.maxE:
                evalueFilter += 1
                continue
            if float(pident) < params.minIdent:
                identFilter +=1
                continue
            if otuId not in hsp:
                hsp[otuId] = {}
                qLen[otuId] = int(qlen)
            if sseqid not in hsp[otuId]:
                hsp[otuId][sseqid] = [(int(qstart), int(qend), float(bitscore))]
                sLen[sseqid] = int(slen)
            else:
                hsp[otuId][sseqid].append((int(qstart), int(qend), float(bitscore)))
        topPerc = params.topPerc/100.0
        with open(output[0], "w") as out, open(log[2], "w") as tLog, open(log[1], "w") as logTax:
            for otuId in hsp.keys():
                hits = []
                for sId, tHsp in hsp[otuId].items():
                    if len(tHsp)>1:
                        used = findTiling(tHsp)
                    else:
                        used = [0]
                    totalLen = 0
                    totalScore = 0
                    for i in used:
                        totalLen += tHsp[i][1] - tHsp[i][0]
                        totalScore += tHsp[i][2]
                    pathStr = ",".join(["%i-%i" % (tHsp[i][0], tHsp[i][1]) for i in used])
                    tLog.write("%s\t%s\t%i\t%i\t%s\t%f\n" % (otuId, sId, len(used), len(tHsp), pathStr, totalScore))
                    if totalLen/min(qLen[otuId], sLen[sId])*100 < params.minCov:
                        covFilter += 1
                        continue
                    linStr = taxDict[sId]
                    hits.append((linStr, totalScore/min(qLen[otuId], sLen[sId]), min(qLen[otuId], sLen[sId])))
                if len(hits) == 0:
                    lineage = "unknown"
                else:
                    sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                    cutoff = 0
                    while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                        cutoff += 1
                    goodHits = [hit[0] for hit in sortedHits[:cutoff]]
                    for h, scr, mlen in sortedHits[:cutoff]:
                        logTax.write("%s\t%s\t%f\t%f\n" % (otuId, h, scr, scr*mlen))
                    lineage = lca(goodHits, params.stringency)
                out.write("%s\t%s\n" % (otuId, lineage))
        with open(log[0], "w") as logOut:
            logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
            logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, params.maxE))
            logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, params.minIdent))
            logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, params.minCov))

rule alignToRdp:
    input: clu="otus/{sample}_{ident}otus_LSU.fasta", db="%(dbFolder)s/rdp_LSU_index.lambda" % config
    output: "lambda/{sample}.{ident}otu_LSU_vs_RDP.m8"
    log: "logs/{sample}_{ident}otuLSU_lambda.log"
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -nm 5000 -p blastn -t {threads} -b -2 -x 40 -as F &> {log}" % config

rule classifyLSU:
    input: lam="lambda/{sample}.{ident}otu_LSU_vs_RDP.m8", tax="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_tax.tsv" % config
    output: "taxonomy/{sample}_{ident}otu_LSU.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/{sample}_LSU_{ident}otuClass.log", "logs/{sample}_{ident}otu_LSU_tax.log", "logs/{sample}_{ident}otu_LSU_tiling.log"
    run:
        taxDict={}
        for line in open(input.tax):
            rId, tax = line.strip().split("\t")
            taxDict[rId] = tax
        classifi = {}
        seqLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
        hsp = {}
        sLen = {}
        qLen = {}
        for line in open(input.lam, encoding="latin-1"):
            total +=1
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.strip().split("\t")
            otuId = qseqid.split("/")[0]
            if float(evalue) > params.maxE:
                evalueFilter += 1
                continue
            if float(pident) < params.minIdent:
                identFilter +=1
                continue
            if otuId not in hsp:
                hsp[otuId] = {}
                qLen[otuId] = int(qlen)
            if sseqid not in hsp[otuId]:
                hsp[otuId][sseqid] = [(int(qstart), int(qend), float(bitscore))]
                sLen[sseqid] = int(slen)
            else:
                hsp[otuId][sseqid].append((int(qstart), int(qend), float(bitscore)))
        topPerc = params.topPerc/100.0
        with open(output[0], "w") as out, open(log[2], "w") as tLog, open(log[1], "w") as logTax:
            for otuId in hsp.keys():
                hits = []
                for sId, tHsp in hsp[otuId].items():
                    if len(tHsp)>1:
                        used = findTiling(tHsp)
                    else:
                        used = [0]
                    totalLen = 0
                    totalScore = 0
                    for i in used:
                        totalLen += tHsp[i][1] - tHsp[i][0]
                        totalScore += tHsp[i][2]
                    pathStr = ",".join(["%i-%i" % (tHsp[i][0], tHsp[i][1]) for i in used])
                    tLog.write("%s\t%s\t%i\t%i\t%s\t%f\n" % (otuId, sId, len(used), len(tHsp), pathStr, totalScore))
                    if totalLen/min(qLen[otuId], sLen[sId])*100 < params.minCov:
                        covFilter += 1
                        continue
                    linStr = taxDict[sId]
                    hits.append((linStr, totalScore/min(qLen[otuId], sLen[sId]), min(qLen[otuId], sLen[sId])))
                if len(hits) == 0:
                    lineage = "unknown"
                else:
                    sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                    cutoff = 0
                    while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                        cutoff += 1
                    goodHits = [hit[0] for hit in sortedHits[:cutoff]]
                    for h, scr, mlen in sortedHits[:cutoff]:
                        logTax.write("%s\t%s\t%f\t%f\n" % (otuId, h, scr, scr*mlen))
                    lineage = lca(goodHits, params.stringency)
                out.write("%s\t%s\n" % (otuId, lineage))
        with open(log[0], "w") as logOut:
            logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
            logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, params.maxE))
            logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, params.minIdent))
            logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, params.minCov))

rule combineCls:
    input: ssu="taxonomy/{sample}_{ident}otu_SSU.class.tsv", its="taxonomy/{sample}_{ident}otu_ITS.class.tsv", lsu="taxonomy/{sample}_{ident}otu_LSU.class.tsv", size="otus/{sample}_{ident}otus.size.tsv"
    output: "taxonomy/{sample}_{ident}_comb.class.tsv"
    run:
        ssu = {}
        for line in open(input.ssu):
            ssuId, ssuTax = line.strip().split("\t")
            ssu[ssuId] = ssuTax
        lsu = {}
        for line in open(input.lsu):
            lsuId, lsuTax = line.strip().split("\t")
            lsu[lsuId] = lsuTax
        its = {}
        for line in open(input.its):
            itsName, itsTax = line.strip().split("\t")
            itsId, itsSizeStr = itsName.strip(";").split(";")
            its[itsId] = itsTax
        with open(output[0], "w") as out:
            for line in open(input.size):
                itsId, sizeStr = line.strip().split("\t")
                size = int(sizeStr)
                itsTax = its.get(itsId, "unknown")
                lsuTax = lsu.get(itsId, "unknown")
                ssuTax = ssu.get(itsId, "unknown")
                out.write("%s\t%i\t%s\t%s\t%s\n" % (itsId, size, ssuTax, itsTax, lsuTax))

rule createOtuTab:
    input: tax="taxonomy/{sampleSet}_{ident}_comb.class.tsv", otu2pClu="otus/{sampleSet}_{ident}otu_preClusterInfo.pic", sample="{sampleSet}_preClu2sample.pic"
    output: "{sampleSet}_otu{ident}_table.tsv"
    run:
        otu2pClu = pickle.load(open(input.otu2pClu, "rb"))
        pCluSample = pickle.load(open(input.sample, "rb"))
        sampleOrder = list(set(pCluSample.values()))
        tab = [["otu", "totalSize", "ssuTax", "itsTax", "lsuTax"] + sampleOrder]
        for line in open(input.tax):
            oId, size, ssuTax, itsTax, lsuTax = line.strip().split("\t")
            pClusters = {}
            for pClu in otu2pClu[oId]:
                tSample = pCluSample[pClu]
                try:
                    pClusters[tSample].append(pClu)
                except KeyError:
                    pClusters[tSample] = [pClu]
            tab.append([oId, size, ssuTax, itsTax, lsuTax])
            for sample in sampleOrder:
                inThisSample = 0
                for pClu in pClusters.get(sample, []):
                    inThisSample += int(pClu.split(";")[1].split("=")[1])
                tab[-1].append(str(inThisSample))
        with open(output[0], "w") as out:
            for line in tab:
                out.write("\t".join(line)+"\n")


rule getCorrectCls:
    input: otuInfo="otus/Lib4-0018_{ident}otuReadInfo.pic", cls="mapping/assignment/Lib4-0018_assignments.tsv"
    output: "taxonomy/Lib4-0018_{ident}otu.mappingClass.tsv"
    run:
        readCls = {}
        for line in open(input.cls):
            read, cls = line.strip().split("\t")
            readCls[read] = cls
        read2otu = pickle.load(open(input.otuInfo, "rb"))
        otuCls = {}
        for read, otu in read2otu.items():
            if otu not in otuCls:
                otuCls[otu] = {}
            tCls = readCls[read]
            try:
                otuCls[otu][tCls] += 1
            except KeyError:
                otuCls[otu][tCls] = 1
        with open(output[0], "w") as out:
            for otu, oCls in otuCls.items():
                out.write("%s\t%s\n" % (otu, ",".join(["%s(%i)" % clsItem for clsItem in oCls.items()])))

rule compareCorrectCls:
    input: correct="taxonomy/Lib4-0018_{ident}otu.mappingClass.tsv", otu= "taxonomy/Lib4-0018_{ident}_comb.class.tsv"
    output: "taxonomy/Lib4-0018_{ident}_combToCorr.class.tsv"
    run:
        corrCls = {}
        for line in open(input.correct):
            otuId, cls = line.strip().split("\t")
            corrCls[otuId] = cls
        with open(output[0], "w") as out:
            for line in open(input.otu):
                otuId, size, ssuTax, itsTax, lsuTax  = line.strip().split("\t")
                corr = corrCls[otuId]
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (otuId, size, corr, ssuTax, itsTax, lsuTax))

rule collectClsStats:
    input: cls="taxonomy/{sampleSet}_97_comb.class.tsv"
    output: complete="taxonomy/{sampleSet}_97_comb.stats.tsv", fungi="taxonomy/{sampleSet}Fungi_97_comb.stats.tsv"
    run:
        tab = []
        for line in open(input.cls):
            oId, size, ssuCls, itsCls, lsuCls = line.strip().split("\t")
            tab.append((oId, int(size), ssuCls, itsCls, lsuCls))

        tab.sort(key=lambda x: x[1], reverse=True)

        with open(output.complete, "w") as out, open(output.fungi, "w") as fungiOut:
            for oId, size, ssuCls, itsCls, lsuCls in tab:
                ssuDepth = 0
                if ssuCls == "unknown":
                    ssuPhyl = "NA"
                else:
                    ssuArr = ssuCls.strip(";").split(";")
                    ssuDepth = len(ssuArr)-1
                    i=0
                    while i<len(ssuArr) and ssuArr[i] != "Fungi":
                        i += 1
                    if i<len(ssuArr):
                        try:
                            ssuPhyl=ssuArr[i+1]
                        except IndexError:
                            ssuPhyl="NA"
                    else:
                        ssuPhyl="non-fungi"
                lsuDepth = 0
                if lsuCls == "unknown":
                    lsuPhyl = "NA"
                else:
                    lsuArr = lsuCls.strip(";").split(";")
                    lsuDepth=len(lsuArr)-1
                    i=0
                    while i<len(lsuArr) and lsuArr[i] != "Fungi":
                        i += 1
                    if i<len(lsuArr):
                        try:
                            lsuPhyl=lsuArr[i+1]
                        except IndexError:
                            lsuPhyl="NA"
                    else:
                        lsuPhyl="non-fungi"
                itsDepth = 0
                if itsCls == "unknown":
                    itsPhyl="NA"
                else:
                    itsArr = itsCls.strip(";").split(";")
                    try:
                        itsPhyl=itsArr[1][3:]
                        itsDepth = len(itsArr)
                    except IndexError:
                        itsPhyl="NA"
                out.write("%s\t%i\tssu\t%s\t%i\n" % (oId, size, ssuPhyl, ssuDepth))
                out.write("%s\t%i\tits\t%s\t%i\n" % (oId, size, itsPhyl, itsDepth))
                out.write("%s\t%i\tlsu\t%s\t%i\n" % (oId, size, lsuPhyl, lsuDepth))
                if ssuDepth>0 and ssuPhyl != "non-fungi":
                    fungiOut.write("%s\t%i\tssu\t%s\t%i\n" % (oId, size, ssuPhyl, ssuDepth))
                    fungiOut.write("%s\t%i\tits\t%s\t%i\n" % (oId, size, itsPhyl, itsDepth))
                    fungiOut.write("%s\t%i\tlsu\t%s\t%i\n" % (oId, size, lsuPhyl, lsuDepth))

rule plotClsComp:
    input: all="taxonomy/{sampleSet}_97_comb.stats.tsv", fungi="taxonomy/{sampleSet}Fungi_97_comb.stats.tsv"
    output: depth="{sampleSet}_clsComp_depth.svg", depthFungi="{sampleSet}_clsComp_depth_fungi.svg", block="{sampleSet}_clsComp_basic.svg"
    run:
        R("""
        library(reshape2)
        library(ggplot2)

        d=read.table("{input.all}", sep="\t")
        colnames(d) = c("oId", "size", "marker", "phylum", "depth")

        d$marker=factor(d$marker, levels=c("ssu", "its", "lsu"))
        d$oId=factor(d$oId, levels=unique(d[order(d$size, decreasing=T),]$oId))
        ggplot(d[1:600,], aes(x=oId,fill=phylum, weight=depth)) + geom_bar() + facet_grid(marker~.) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
        ggsave("{output.depth}", width=16, height=10)


        f=read.table("{input.fungi}", sep="\t")
        colnames(f) = c("oId", "size", "marker", "phylum", "depth")
        f$marker=factor(f$marker, levels=c("ssu", "its", "lsu"))
        
        f$oId=factor(f$oId, levels=unique(f[order(f$marker, f$phylum, f$depth),]$oId))
        ggplot(f[1:600,], aes(x=oId,fill=phylum, weight=depth)) + geom_bar() + facet_grid(marker~.) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6))
        ggsave("{output.depthFungi}", width=16, height=10)


        d$marker=factor(d$marker, levels=c("lsu", "ssu", "its"))
        ggplot(d[1:600,], aes(x=oId, y=marker, fill=phylum)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
        ggsave("{output.block}", width=16, height=10)
        """)

rule plotDephtOverview:
    input: "taxonomy/all_97_comb.stats.tsv"
    output: "classificationDepth.svg"
    params: minSize=5
    run:
        R("""
        library(ggplot2)
        
        d=read.table("{input}", sep="\t")
        colnames(d) = c("oId", "size", "marker", "phylum", "depth")
        d$marker=factor(d$marker, levels=c("ssu", "its", "lsu"))
        
        ggplot(d[d$size>={params.minSize},]) + geom_bar(aes(factor(depth), weight=size)) + facet_grid(.~marker) + scale_x_discrete(labels=c("NA", "kingdom", "phylum", "class", "order", "family", "genus", "species")) + coord_flip() + labs(y="read count", x="classfied to") + theme_bw()
        ggsave("{output}", width=16, height=10)
        """)

rule compareCls:
    input: cls="taxonomy/all_97_comb.class.tsv"
    output: diff="all_clsDiff.tsv", comp="all_clsComp.tsv"
    run:
        ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        diff = {"ssu_lsu": {}, "ssu_its": {}, "its_lsu":{}}
        size = {}
        allCls = {}
        com = {}
        for line in open(input.cls):
            oId, tSize, ssuCls, itsCls, lsuCls = line.strip().split("\t")
            cls = {}
            maxLen=0
            if ssuCls == "unknown":
                cls["ssu"] = None
            else:
                cls["ssu"] = ssuCls.split(";")
                if len(cls["ssu"]) > maxLen:
                    maxLen=len(cls["ssu"])
            if lsuCls == "unknown":
                cls["lsu"] = None
            else:
                cls["lsu"] = lsuCls.split(";")
                if len(cls["lsu"]) > maxLen:
                    maxLen=len(cls["lsu"])
            if itsCls == "unknown":
                cls["its"] = None
            else:
                cls["its"] = ["Eukaryota"] + [c.split("__", 1)[-1] for c in itsCls.split(";")]
                if len(cls["its"]) > maxLen:
                    maxLen=len(cls["its"])
            allCls[oId] = cls
            com[oId] = []
            
            for lvl in range(maxLen):
                tCls = [None, None, None]
                if not cls["ssu"] is None and len(cls["ssu"]) > lvl:
                    tCls[0] = cls["ssu"][lvl]
                if not cls["its"] is None and len(cls["its"]) > lvl:
                    tCls[1] = cls["its"][lvl]
                if not cls["lsu"] is None and len(cls["lsu"]) > lvl:
                    tCls[2] = cls["lsu"][lvl]
                if not tCls[0] is None and (tCls[0] == tCls[1] or tCls[0] == tCls[2]):
                    com[oId].append(tCls[0])
                elif not tCls[1] is None and tCls[1] == tCls[2]:
                    com[oId].append(tCls[1])
                else:
                    com[oId].append(None)
                
            for mrk1, mrk2 in [("ssu","lsu"), ("ssu", "its"), ("its", "lsu")]:
                if cls[mrk1] is None or cls[mrk2] is None:
                    continue
                for r, rank in enumerate(ranks):
                    if r >= len(cls[mrk1]) or r >= len(cls[mrk2]):
                        break
                    if cls[mrk1][r] != cls[mrk2][r]:
                        if rank not in diff["%s_%s" % (mrk1, mrk2)]:
                            diff["%s_%s" % (mrk1, mrk2)][rank] = {}
                        try:
                            diff["%s_%s" % (mrk1, mrk2)][rank]["%s<->%s" % (cls[mrk1][r], cls[mrk2][r])].append(oId)
                        except KeyError:
                            diff["%s_%s" % (mrk1, mrk2)][rank]["%s<->%s" % (cls[mrk1][r], cls[mrk2][r])] = [oId]

        with open(output.diff, "w") as out:
            for comp in diff.keys():
                mrk1, mrk2 = comp.split("_")
                for rank, diffData in diff[comp].items():
                    for entry, otus in diffData.items():
                        out.write("%s\t%s\t%s\t%s\t%s\t%i\n" % (comp, mrk1, mrk2, rank, entry, len(otus)))
        with open(output.comp, "w") as out:
            for oId in com.keys():
                out.write("%s\t%s;" % (oId, ";".join([str(c) for c in com[oId]])))
                for marker in ["ssu", "its", "lsu"]:
                    out.write("\t")
                    cls = allCls[oId][marker]
                    if cls is None:
                        out.write("NA;")
                    else:
                        for r in range(len(cls)):
                            if r < len(com[oId]):
                                if cls[r] != com[oId][r]:
                                    out.write("%s;" % cls[r])
                                else:
                                    out.write("--;")
                            else:
                                out.write("%s;" % cls[r])
                out.write("\n")
                

rule combineIsolateCls:
    input: expand("taxonomy/{spec}_97_comb.class.tsv", spec=isolates.keys())
    output: "taxonomy/isolates_97_comb.class.tsv"
    run:
        with open(output[0], "w") as out:
            for inputFile in input:
                spec = inputFile.split("/", 1)[1].split("_", 1)[0]
                for line in open(inputFile):
                    out.write("%s\t%s" % (spec, line))
            

def findTiling(hsp):
    G=nx.DiGraph()
    nodes = []
    b_edges = []
    for h, tHsp in enumerate(hsp):
        nodes.extend(["S%i" % h, "E%i" % h])
        if tHsp[0] < tHsp[1]:
            #forward
            b_edges.append(("S%i" % h, "E%i" % h, -tHsp[2])) #use negative bit score as edge weight to be able to use minimum path algorithm to find maximum path
        else:
            #reverse
            b_edges.append(("E%i" % h, "S%i" % h, -tHsp[2]))
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(b_edges)
    #create edges between all non-overlapping matches that are in the same direction
    c_edges = []
    for i in range(len(hsp)):
        for j in range(len(hsp)):
            if i==j:
                continue
            if (hsp[i][0] < hsp[i][1]) and (hsp[j][0] < hsp[j][1]) and (hsp[i][1] < hsp[j][0]):
                #forward compatibility edge
                c_edges.append(("E%i" % i, "S%i" % j, 0))
            elif (hsp[i][0] > hsp[i][1]) and (hsp[j][0] > hsp[j][1]) and (hsp[i][1] > hsp[j][0]):
                #reverse compatibility edge
                c_edges.append(("S%i" % j, "E%i" % i, 0))
    G.add_weighted_edges_from(c_edges)

    #find sources and start node
    s_edges = []
    for node, id in G.in_degree_iter():
        if id==0:
            s_edges.append(("START", node, 0))
    G.add_weighted_edges_from(s_edges)
    #find sinks and add end node
    e_edges = []
    for node, od in G.out_degree_iter():
        if od==0:
            e_edges.append((node, "END", 0))
    G.add_weighted_edges_from(e_edges)
    #find shortest (ie. maximum score) path from start to end
    pre, dist = nx.floyd_warshall_predecessor_and_distance(G)
    if dist["END"]["START"] < dist["START"]["END"]:
        start = "END"
        end = "START"
    else:
        start = "START"
        end = "END"
    current = end
    path = []
    try:
        while pre["START"][current] != start:
            path.append(pre["START"][current])
            current = pre["START"][current]
    except KeyError:
        print(hsp)
        print(G.edges())
        print(path)
        print(pre)
        raise
    
    return [int(p[1:]) for p in path[::-2]]

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


#####################################################################

rule indiDerep:
    input: "itsx/{sample}.{marker}.fasta"
    output: fasta="indiDerep/{sample}_{marker}.derep.fasta", info="indiDerep/{sample}_{marker}.repseq.pic", txt="indiDerep/{sample}_{marker}.uc.txt"
    log: "logs/{sample}_{marker}_derep.log"
    run:
        shell("%(vsearch)s --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizeout --log {log} &> /dev/null" % config)
        repseq = {}
        for line in open(output.txt):
            arr = line.strip().split("\t")
            if arr[0] == "C":
                pass
            elif arr[0] == "S":
                seq = arr[-2]
                if seq in repseq:
                    raise ValueError("Starting existing cluster")
                repseq[seq] = [seq]
            elif arr[0] == "H":
                seq, cluster = arr[-2:]
                repseq[cluster].append(seq)
            else:
                raise ValueError("Unknown record type: %s" % arr[0])
        pickle.dump(repseq, open(output.info, "wb"))

rule indiCluster:
    input: "indiDerep/{sample}_{marker}.derep.fasta"
    output: fasta="indiCluster/{sample}_{marker}_clu.fasta", uc="indiCluster/{sample}_{marker}_clu.uc.tsv"
    log: "logs/{sample}_indiClustering_{marker}.log"
    threads: 6
    shell:
        "%(vsearch)s --cluster_size {input} --id 0.97 --sizein --sizeout --centroids {output.fasta} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

rule cp58s:
    input: "itsx/{sample}.5_8S.fasta"
    output: "itsx/{sample}.58S.fasta"
    shell:
        "cp {input} {output}"

rule indiCluReads:
    input: clsInfo="indiCluster/{sample}_{marker}_clu.uc.tsv", repseq="indiDerep/{sample}_{marker}.repseq.pic"
    output: info="indiCluster/{sample}_{marker}_cluInfo.pic"
    run:
        repseq=pickle.load(open(input.repseq, "rb"))
        clusterReads={}
        tReads = None
        repSeq = None
        for line in open(input.clsInfo):
            arr = line.strip().split("\t")
            if arr[0] == "C":
                pass
            elif arr[0] == "S":
                seq = arr[-2]
                if seq in repseq:
                    raise ValueError("Starting existing cluster")
                seqId = seq.split(";", 1)[0]
                clusterReads[seqId] = [s.split("|", 1)[0] for s in repseq[seqId]]
            elif arr[0] == "H":
                seq, cluster = arr[-2:]
                seqId = seq.split(";", 1)[0]
                clusterId = cluster.split(";", 1)[0]
                clusterReads[clusterId].extend([s.split("|", 1)[0] for s in repseq[seqId]])
            else:
                raise ValueError("Unknown record type: %s" % arr[0])
        with open(output.info, "wb") as pOut:
            pickle.dump(clusterReads, pOut)

rule indiCluOverlap:
    input: ssu="indiCluster/{sample}_SSU_cluInfo.pic", its1="indiCluster/{sample}_ITS1_cluInfo.pic", r58s="indiCluster/{sample}_58S_cluInfo.pic", its2="indiCluster/{sample}_ITS2_cluInfo.pic", lsu="indiCluster/{sample}_LSU_cluInfo.pic"
    output: nodes="clusterGraph/{sample}_clusterGraphNodes.tsv", edges="clusterGraph/{sample}_clusterGraphEdges.tsv"
    params: minCluSize=2, minCluOverlap=2
    run:
        ssu=pickle.load(open(input.ssu, "rb"))
        its1=pickle.load(open(input.its1, "rb"))
        r58s=pickle.load(open(input.r58s, "rb"))
        its2=pickle.load(open(input.its2, "rb"))
        lsu=pickle.load(open(input.lsu, "rb"))
        nodes = {}
        for markerData in [ssu, its1, r58s, its2, lsu]:
                for cId, reads in markerData.items():
                    if len(reads) >= params.minCluSize:
                        nodes[cId] = reads
        with open(output.nodes, "w") as out:
            out.write("cluster\ttype\tsize\n")
            for cId, reads in nodes.items():
                out.write("%s\t%s\t%i\n" % (cId, cId.split("|")[2], len(reads)))
                    
        edges = []
        for start, end in [(ssu, its1), (its1, r58s), (r58s, its2), (its2, lsu)]:
            for startId in start:
                if not startId in nodes:
                    continue
                for endId in end:
                    if not endId in nodes:
                        continue
                    weight = len(set(start[startId]) & set(end[endId]))
                    if weight >= params.minCluOverlap:
                        edges.append((startId, endId, weight))
        with open(output.edges, "w") as out:
            out.write("start\tend\tweight\n")
            for e in edges:
                out.write("%s\t%s\t%i\n" % e)

rule findComponents:
    input: edges="clusterGraph/{sample}_clusterGraphEdges.tsv"
    output: comp="clusterGraph/{sample}_components.txt"
    run:
        G=nx.Graph()
        with open(input.edges) as inStream:
            _ = next(inStream) # header
            for line in inStream:
                start, end, weight = line.strip().split("\t")
                G.add_edge(start, end)
            with open(output.comp, "w") as out:
                for comp in nx.connected_components(G):
                    out.write("\t".join(comp)+"\n")

rule indiCluClass:
    input: ssu="indiCluster/{sample}_SSU_cluInfo.pic", its1="indiCluster/{sample}_ITS1_cluInfo.pic", r58s="indiCluster/{sample}_58S_cluInfo.pic", its2="indiCluster/{sample}_ITS2_cluInfo.pic", lsu="indiCluster/{sample}_LSU_cluInfo.pic", cls="mapping/assignment/{sample}_assignments.tsv", comp="clusterGraph/{sample}_components.txt"
    output: tab="clusterGraph/{sample}_components_class.tsv", lab="clusterGraph/{sample}_clusterGraphCls.tsv"
    run:
        reads = {
        "SSU": pickle.load(open(input.ssu, "rb")),
        "ITS1": pickle.load(open(input.its1, "rb")),
        "5.8S": pickle.load(open(input.r58s, "rb")),
        "ITS2": pickle.load(open(input.its2, "rb")),
        "LSU": pickle.load(open(input.lsu, "rb"))
        }
        readCls = {}
        for line in open(input.cls):
            read, cls = line.strip().split("\t")
            readCls[read] = cls
        with open(output.tab, "w") as out, open(output.lab, "w") as labOut:
            labOut.write("nodeName\tclassification\n")
            for cNr, line in enumerate(open(input.comp)):
                comp = line.strip().split("\t")
                for clu in comp:
                    sId, _, marker = clu.rsplit("|")
                    cluCls = {}
                    for read in reads[marker][clu]:
                        tCls = readCls[read]
                        try:
                            cluCls[tCls] += 1
                        except KeyError:
                            cluCls[tCls] = 1
                    for cls, count in cluCls.items():
                        out.write("Comp%i\t%s\t%s\t%s\t%i\n" % (cNr, marker, clu, cls, count))
                    labOut.write("%s\t%s\n" % (clu, "+".join(["%s(%i)" % i for i in cluCls.items()])))

rule componentAlign:
    input: comp="clusterGraph/{sample}_components.txt", fasta="itsx/{sample}.{marker}.fasta", info="indiCluster/{sample}_{marker}_cluInfo.pic"
    output: dynamic("clusterAln/{sample}_compAln_{comp}_{marker}.fasta")
    log: comp="logs/{sample}_components.log", aln="logs/{sample}_cluAlign.log"
    threads: 6
    run:
        with open(log.comp, "w") as logFile:
            seqRecs = {}
            for rec in SeqIO.parse(open(input.fasta), "fasta"):
                seqRecs[rec.id.split("|", 1)[0]] = rec
            readInfo=pickle.load(open(input.info, "rb"))
            for l,line in enumerate(open(input.comp)):
                comp = line.strip().split("\t")
                markers = {}
                for clu in comp:
                    try:
                        markers[clu.rsplit("|", 1)[-1]].append(clu)
                    except KeyError:
                        markers[clu.rsplit("|", 1)[-1]] = [clu]
                
                baseFileName = "clusterAln/tmp_comp%i_%s" % (l, wildcards.marker)
                #write all reads in the component to temp file
                readFilePath = "%s_reads.fasta" % baseFileName
                alignmentPath="clusterAln/%s_compAln_%i_%s.fasta" % (wildcards.sample, l, wildcards.marker)
                if wildcards.marker=="58S":
                    marker="5.8S"
                else:
                    marker=wildcards.marker
                with open(readFilePath, "w") as readFile:
                    if marker not in markers:
                        logFile.write("No %s for component %i" % (marker, l))
                        #create dummy file 
                        open(alignmentPath, "w")
                        continue
                    for c, clu in enumerate(markers[marker]):
                        for readId in readInfo[clu]:
                            rec = seqRecs[readId]
                            rec.id = "%s|%s" % (c, rec.id)
                            readFile.write(rec.format("fasta"))
                #run alignment
                shell("%(mafft)s" % config + " --op 0.5 --ep 0.5 --thread {threads} %s > %s 2> {log.aln}" % (readFilePath, alignmentPath))
                shell("rm clusterAln/tmp_*")

def fastaConsensus(inFasta, warn=0.75, log=sys.stderr):
    prof = []
    for rec in SeqIO.parse(inFasta, "fasta"):
        for c, char in enumerate(rec.seq):
            try:
                prof[c][char.upper()] += 1
            except IndexError:
                prof.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
                prof[c][char.upper()] += 1
    cons = []
    for i, col in enumerate(prof):
        cov = float(sum(col.values()))
        most = sorted(list(col.items()), key=lambda x: x[1], reverse=True)[0]
        if most[1]/cov < warn:
            log.write("At position %i the most common base %s only has a frequency of %f\n" % (i, most[0], most[1]/cov))
        if most[0] != "-":
            cons.append(most[0])
    return "".join(cons)

#def samConsensus(inSam, warn=0.75, log=sys.stderr):
#    prof = []
#    for line in inSam:
#        if line[0] == "@":
#            continue
#        aSeq = []
#        readId, flag, refId, pos, mapQ, cigar, pnext, rnext, alen, seq, seqQ, tags = line.strip().split("\t", 11)
#        aSeq.append("-"*(int(pos)-1))
#        cur=0
#        for entry in re.findall("[0-9]+[MIDNSHP=X]", cigar):
#            l = int(entry[:-1])
#            if entry[-1] == "M":
#                aSeq.append(seq[cur:cur+l])
#                cur += l
#            elif entry[-1] == "I":
#                pass
#                cur += l
#            elif entry[-1] == "D":
#                aSeq.append("-"*l)
#            else:
#                raise NotImplementedError("Cigar entry %s not implemented yet" % entry[-1])
##        log.write(">%s\n%s\n" % (readId, "".join(aSeq)))
#        for c,char in enumerate("".join(aSeq)):
#            try:
#                prof[c][char] += 1
#            except IndexError:
#                prof.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
#                prof[c][char] += 1
#    cons = []
#    for i, col in enumerate(prof):
#        cov = float(sum(col.values()))
#        most = sorted(list(col.items()), key=lambda x: x[1], reverse=True)[0]
#        if most[1]/cov < warn:
#            log.write("At position %i the most common base %s only has a frequency of %f\n" % (i, most[0], most[1]/cov))
#        if most[0] != "-":
#            cons.append(most[0])
#    return "".join(cons)


