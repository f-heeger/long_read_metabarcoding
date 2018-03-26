import pickle
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from snakemake.utils import min_version, R

import networkx as nx

shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow samba

configfile: "config.json"

sampleInfo = json.load(open("samples.json"))

samples = [sId for sId, sType in sampleInfo["sampleType"].items() if sType == "env"]
isolates = {}
for sId, sType in sampleInfo["sampleType"].items():
    if sType == "isolate":
        try:
            isolates[sampleInfo["sampleName"][sId]].append(sId)
        except KeyError:
            isolates[sampleInfo["sampleName"][sId]] = [sId]

allSamples = list(sampleInfo["sampleName"].keys())

sampleName = sampleInfo["sampleName"]

mockSamples = [sId for sId, sType in sampleInfo["sampleType"].items() if sType == "mock"]

####################################################################
# includes

include: "createDBs.snakemake.py"
include: "readProcessing.snakemake.py"
include: "mapping.snakemake.py"
include: "chimera_analysis.snakemake.py"

####################################################################

rule all:
    input: "taxonomy/all_97_comb.class.tsv", "all_clsComp_depth.svg", "all_clsComp_depth_fungi.svg", "all_clsComp_basic.svg", "taxonomy/Lib4-0018_97_combToCorr.class.tsv", "chimeraCyclesRelativeBarplot.svg", "chimera_comp_sankey.svg", expand(["mapping/{stage}MockComp.svg", "mapping/{stage}ErrorRates.svg"], stage=["raw", "filtered"]), "readNumbers.svg"

rule concatItsxResult:
    """Concatenate ITSx results from different samples"""
    input: expand("itsx/{sample}.{{marker}}.fasta", sample=samples)
    output: "catItsx/all.{marker}.fasta"
    shell:
        "cat {input} > {output}"

def catIsolatesInput(wildcards):
    """determine input data for catIsolate rule"""
    return ["itsx/%s.%s.fasta" % (s, wildcards.marker) for s in isolates[wildcards.spec]]

rule concatIsolates:
    """Concatenate ITSx results for an isolate sample"""
    input: catIsolatesInput
    output:"isolates/{spec}.{marker}.fasta"
    shell:
        "cat {input} > {output}"

def otuInput(wildcards):
    """determine input file for OTU clustering according to sample wildcard"""
    if wildcards.sample == "all":
        return "catItsx/all.full.fasta"
    elif wildcards.sample in isolates:
        return "isolates/%s.full.fasta" % wildcards.sample
    else:
        return "itsx/%s.full.fasta" % wildcards.sample

rule otuCluster:
    """Cluster OTUs by ITS sequence with vsearch"""
    input: otuInput
    output: fasta="otus/{sample}_{ident}otus.fasta", uc="otus/{sample}_{ident}otus.uc.tsv"
    log: "logs/{sample}_{ident}otuClustering.log"
    threads: 3
    shell:
        "%(vsearch)s --cluster_size {input} --relabel otu --sizein --sizeout --iddef 0 --id 0.{wildcards.ident} --minsl 0.9 --centroids {output.fasta} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

def otuReadsInput(wildcards):
    """determine input data for otuReads rule according to sample wildcard"""
    if wildcards.sample == "all":
        return ["otus/all_%sotus.uc.tsv" % wildcards.ident] + expand("preclusters/{sample}_cluInfo.tsv", sample=samples)
    elif wildcards.sample in isolates:
        return ["otus/%s_%sotus.uc.tsv" % (wildcards.sample, wildcards.ident)] + expand("preclusters/{sample}_cluInfo.tsv", sample=isolates[wildcards.sample])
    else:
        return  ["otus/%s_%sotus.uc.tsv" % (wildcards.sample, wildcards.ident), "preclusters/%s_cluInfo.tsv" % (wildcards.sample)]

rule otuReads:
    """Generate different infos for the OTUs (as dictonaries in pickle files)
    
    *otus.size.tsv: number of reads in the OTU
    *otuReadInfo.pic: reads for each OTU
    *otu_preClusterInfo.pic: preClusters for each OTU
    *otu_repSeq.pic: representative sequence of the OTU
    """
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
    """determine input data for transferOtus rule according to sample wildcard"""
    if wildcards.sample == "all":
        return ["otus/all_%sotu_repSeq.pic" % wildcards.ident, "catItsx/all.%s.fasta" % wildcards.marker]
    elif wildcards.sample in isolates:
        return ["otus/%s_%sotu_repSeq.pic" % (wildcards.sample, wildcards.ident), "isolates/%s.%s.fasta" % (wildcards.sample, wildcards.marker)]
    else:
        return ["otus/%s_%sotu_repSeq.pic" % (wildcards.sample, wildcards.ident), "itsx/%s.%s.fasta" % (wildcards.sample, wildcards.marker)]


rule transferOtus:
    """Create file with non-ITS marker sequence (SSU or LSU) per OTU"""
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
    """Search for local sequence matches of ITS sequences in the UNITE database with lambda"""
    input: clu="otus/{sample}_{ident}otus.fasta", db="%(dbFolder)s/UNITE_%(uniteVersion)s.index.lambda" % config
    output: "lambda/{sample}.{ident}otu_vs_UNITE.m8"
    log: "logs/{sample}_{ident}otu_lambda.log"
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -p blastn -t {threads} &> {log}" % config

rule classifyITS:
    """Classify OTU based on ITS matches with LCA approach"""
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
    """Search for local sequence matches of SSU sequences in the SILVA database with lambda"""
    input: clu="otus/{sample}_{ident}otus_SSU.fasta", db="%(dbFolder)s/silva_SSU_index.lambda" % config
    output: "lambda/{sample}.{ident}otu_SSU_vs_SILVA.m8"
    log: "logs/{sample}_{ident}otu_SSU_lambda.log"
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -nm 20000 -p blastn -t {threads} -b -2 -x 30 -as F &> {log}" % config

rule classifySSU:
    """Classify OTU based on SSU matches with LCA approach"""
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
                    tLog.write("%s\t%s\t%i\t%i\t%s\t%f\t%f\n" % (otuId, sId, len(used), len(tHsp), pathStr, totalScore, totalScore/min(qLen[otuId], sLen[sId])))
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
    """Search for local sequence matches of LSU sequences in the RDP LSU database with lambda"""
    input: clu="otus/{sample}_{ident}otus_LSU.fasta", db="%(dbFolder)s/rdp_LSU_index.lambda" % config
    output: "lambda/{sample}.{ident}otu_LSU_vs_RDP.m8"
    log: "logs/{sample}_{ident}otuLSU_lambda.log"
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -nm 5000 -p blastn -t {threads} -b -2 -x 40 -as F &> {log}" % config

rule classifyLSU:
    """Classify OTU based on LSU matches with LCA approach"""
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
                    tLog.write("%s\t%s\t%i\t%i\t%s\t%f\t%f\n" % (otuId, sId, len(used), len(tHsp), pathStr, totalScore, totalScore/min(qLen[otuId], sLen[sId])))
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
    """Create table with OTU classifications with SSU, ITS, LSU"""
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
    """Create a OTU table, with size, classiifcation and occurence information"""
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
    """Get "correct" classification for OTUs based on reference mapping"""
    input: otuInfo="otus/Lib4-0018_{ident}otuReadInfo.pic", cls="mapping/assignment/Lib4-0018_filtered_assignments.tsv"
    output: "taxonomy/Lib4-0018_{ident}otu.mappingClass.tsv"
    run:
        readCls = {}
        for line in open(input.cls):
            read, cls = line.strip().split("\t")[:2]
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
    """Create table with OTU classifications including the "correct" one"""
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
    """Create table of classifications (including depth) as input for ggplot"""
    input: cls="taxonomy/{sampleSet}_97_comb.class.tsv"
    output: complete="taxonomy/{sampleSet}_97_comb.stats.tsv" 
    run:
        ranks = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        tab = []
        for line in open(input.cls):
            oId, size, ssuCls, itsCls, lsuCls = line.strip().split("\t")
            tab.append((oId, int(size), ssuCls, itsCls, lsuCls))

        tab.sort(key=lambda x: x[1], reverse=True)

        with open(output.complete, "w") as out:
            for oId, size, ssuClsStr, itsClsStr, lsuClsStr in tab:
                ssuDepth = 0
                ssuCls = dict(zip(ranks, ["NA"]*8))
                if ssuClsStr != "unknown":
                    ssuArr = ssuClsStr.strip(";").split(";")
                    ssuDepth = len(ssuArr)-1
                    for r, rank in enumerate(ranks):
                        try:
                            ssuCls[rank] = ssuArr[r]
                        except IndexError:
                            break
                lsuDepth = 0
                lsuCls = dict(zip(ranks, ["NA"]*8))
                if lsuClsStr != "unknown":
                    lsuArr = lsuClsStr.strip(";").split(";")
                    lsuDepth = len(lsuArr)-1
                    for r, rank in enumerate(ranks):
                        try:
                            lsuCls[rank] = lsuArr[r]
                        except IndexError:
                            break
                itsDepth = 0
                itsCls = dict(zip(ranks, ["NA"]*8))
                if itsClsStr != "unknown":
                    itsArr = itsClsStr.strip(";").split(";")
                    itsDepth = len(itsArr)
                    itsCls["kingdom"] = "Eukaryota"
                    for r, rank in enumerate(ranks[1:]):
                        try:
                            itsCls[rank] = itsArr[r][3:]
                        except IndexError:
                            break
                out.write("%s\t%i\tssu\t%s\t%i\n" % (oId, size, "\t".join([ssuCls[r] for r in ranks]), ssuDepth))
                out.write("%s\t%i\tits\t%s\t%i\n" % (oId, size, "\t".join([itsCls[r] for r in ranks]), itsDepth))
                out.write("%s\t%i\tlsu\t%s\t%i\n" % (oId, size, "\t".join([lsuCls[r] for r in ranks]), lsuDepth))

rule clsSummary:
    """Collect some summary stats of how many OTUs were assigend to taxonomic 
    ranks for the paper abstract"""
    input: "taxonomy/{sampleSet}_97_comb.stats.tsv"
    output: "taxonomy/{sampleSet}_97_clsStats.tsv"
    run:
        R("""
        d=read.table("{input}", sep="\t")
        colnames(d) = c("oId", "size", "mrk", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "depth")
        a=aggregate(depth ~ oId, subset(d, size>1), max)
        ranks=c("kingdom", "phylum", "class", "order", "family", "genus", "species")
        s = data.frame(rank=numeric(0), nubmer=numeric(0))
        for (i in 1:7) {{
            s = rbind(s, data.frame(rank=ranks[i], number=sum(a$depth>=i)))
        }}
        write.table(s, "{output}", sep="\t", row.names=F)
        """)

rule plotClsComp:
    """Create plots of classifications depth"""
    input: all="taxonomy/{sampleSet}_97_comb.stats.tsv"
    output: depth="{sampleSet}_clsComp_depth.svg", depthFungi="{sampleSet}_clsComp_depth_fungi.svg", block="{sampleSet}_clsComp_basic.svg"
    run:
        R("""
        library(reshape2)
        library(ggplot2)

        d=read.table("{input.all}", sep="\t")
        colnames(d) = c("oId", "size", "marker", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "depth")

        d$marker=factor(d$marker, levels=c("ssu", "its", "lsu"))
        d$oId=factor(d$oId, levels=unique(d[order(d$size, decreasing=T),]$oId))
        ggplot(d[1:600,], aes(x=oId,fill=phylum, weight=depth)) + geom_bar() + facet_grid(marker~.) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
        ggsave("{output.depth}", width=16, height=10)

        fo = unique(d[d$kingdom=="Fungi",]$oId)
        f=subset(d, d$oId %in% fo)
        
        f$oId=factor(f$oId, levels=unique(f[order(f$marker, f$phylum, -f$depth),]$oId))
         ggplot(f[1:600,], aes(x=oId,fill=phylum, weight=depth)) + geom_bar() + facet_grid(marker~.) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6)) + scale_y_continuous(breaks=seq(1, 7, 1), labels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
        ggsave("{output.depthFungi}", width=16, height=10)

        #sum(table(subset(f[1:600,], marker=="ssu" & size>1)$depth)[c("2","3","4","5","6","7")])
        #sum(table(subset(f[1:600,], marker=="lsu" & size>1)$depth)[c("2","3","4","5","6","7")])
        #sum(table(subset(f[1:600,], marker=="its" & size>1)$depth)[c("2","3","4","5","6","7")])

        d$marker=factor(d$marker, levels=c("lsu", "ssu", "its"))
        ggplot(d[1:600,], aes(x=oId, y=marker, fill=phylum)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
        ggsave("{output.block}", width=16, height=10)
        """)


rule compareCls:
    """Compare classifications from different markers for plotting with ggplot
    
    Is a OTU that is classified by marker A at a certain level also classified 
    by marker B and C at this level.
    """
    input: cls="taxonomy/{sampleSet}_97_comb.class.tsv"
    output: diff="{sampleSet}_clsDiff.tsv", comp="{sampleSet}_clsComp.tsv", diffStat="{sampleSet}_clsDiffStat.tsv"
    run:
        ranks=["eukaryota", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        diff = {"ssu_lsu": {}, "ssu_its": {}, "its_lsu":{}}
        size = {}
        allCls = {}
        com = {}
        stat = {}
        fungi = {}
        for line in open(input.cls):
            oId, tSize, ssuCls, itsCls, lsuCls = line.strip().split("\t")
            size[oId] = int(tSize)
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
                cls["lsu"] = [lsuEntry.replace(" ", "_") for lsuEntry in lsuCls.split(";")]
                if len(cls["lsu"]) > maxLen:
                    maxLen=len(cls["lsu"])
            if itsCls == "unknown":
                cls["its"] = None
            else:
                cls["its"] = ["Eukaryota"] + [c.split("__", 1)[-1] for c in itsCls.split(";")]
                if len(cls["its"]) > maxLen:
                    maxLen=len(cls["its"])
            fungi[oId] = any([cls["ssu"] and len(cls["ssu"]) > 1 and cls["ssu"][1] == "Fungi", 
                              cls["its"] and len(cls["its"]) > 1 and cls["its"][1] == "Fungi", 
                              cls["lsu"] and len(cls["lsu"]) > 1 and cls["lsu"][1] == "Fungi"])
            allCls[oId] = cls
            com[oId] = []
            stat[oId] = {}
            if not cls["ssu"] is None:
                stat[oId]["ssu"] = {"its": ["NA"]*8,
                                    "lsu": ["NA"]*8}
            if not cls["its"] is None:
                stat[oId]["its"] = {"ssu": ["NA"]*8,
                                    "lsu": ["NA"]*8}
            if not cls["lsu"] is None:
                stat[oId]["lsu"] = {"ssu": ["NA"]*8,
                                    "its": ["NA"]*8}
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
                for i, mrk in [(0, "ssu"), (1, "its"), (2, "lsu")]:
                    if not tCls[i] is None:
                        (o1, oMrk1), (o2, oMrk2) = set([(0, "ssu"), (1, "its"), (2, "lsu")]) - set([(i,mrk)])
                        if tCls[o1] is None:
                            stat[oId][mrk][oMrk1][lvl] = "blank"
                        else:
                            if tCls[i] == tCls[o1]:
                                stat[oId][mrk][oMrk1][lvl] = "same"
                            else:
                                stat[oId][mrk][oMrk1][lvl] = "diff"
                        if tCls[o2] is None:
                            stat[oId][mrk][oMrk2][lvl] = "blank"
                        else:
                            if tCls[i] == tCls[o2]:
                                stat[oId][mrk][oMrk2][lvl] = "same"
                            else:
                                stat[oId][mrk][oMrk2][lvl] = "diff"
                    
            for mrk1, mrk2 in [("ssu","lsu"), ("ssu", "its"), ("its", "lsu")]:
                if cls[mrk1] is None or cls[mrk2] is None:
                    continue
                for r, rank in enumerate(ranks):
                    if r >= len(cls[mrk1]) or r >= len(cls[mrk2]):
                        break
                    if cls[mrk1][r] != cls[mrk2][r]:
                        if r == 0:
                            parent="root"
                        else:
                            parent=cls[mrk1][r-1]
                        if r>2:
                            pphylum=cls[mrk1][2]
                        else:
                            pphylum="NA"
                        if rank not in diff["%s_%s" % (mrk1, mrk2)]:
                            diff["%s_%s" % (mrk1, mrk2)][rank] = []
                        diff["%s_%s" % (mrk1, mrk2)][rank].append((oId, cls[mrk1][r], cls[mrk2][r], parent, pphylum))
                        break

        with open(output.diff, "w") as out:
            for comp in diff.keys():
                mrk1, mrk2 = comp.split("_")
                for rank, diffData in diff[comp].items():
                    for oId, m1Cls, m2Cls, parent, phylum in diffData:
                        out.write("%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\n" % (comp, mrk1, mrk2, rank, oId, size[oId], m1Cls, m2Cls, parent, phylum))
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
                                elif not com[oId][r] is None:
                                    out.write("--;")
                            else:
                                out.write("%s;" % cls[r])
                out.write("\n")
        with open(output.diffStat, "w") as sOut:
            for oId, data in stat.items():
                for mrk, mrkData in data.items():
                    for oMrk, values in mrkData.items():
                        sOut.write("%s\t%i\t%s\t%s\t%s\t%s\n" % (oId, size[oId], fungi[oId], mrk, oMrk, "\t".join(values)))

rule plotDiff:
    """Plot comparison data of classification with different markers"""
    input: "{sampleSet}_clsDiffStat.tsv"
    output: bars="{sampleSet}_clsDiffStat.svg", prop="{sampleSet}_clsDiffPropSame.svg"
    run:
        R("""
        library(ggplot2)
        library(reshape2)

        d=read.table("{input}", sep="\t")
        colnames(d) = c("oId", "size", "fungi", "mrk1", "mrk2", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

        for (lvl in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {{
            d[,lvl] = factor(d[,lvl], levels=c("same", "diff", "blank", NA))
        }}

        m=melt(d, id.vars=c("oId", "size", "fungi", "mrk1", "mrk2"))
        m$mrk1=factor(m$mrk1, levels=c("ssu", "its", "lsu"))
        m$mrk2=factor(m$mrk2, levels=c("ssu", "its", "lsu"))
        m$value = factor(m$value, levels=c(NA, "blank", "diff", "same"))
        
        m$mrk1=factor(m$mrk1, levels=c("ssu", "its", "lsu"))
        m$mrk2=factor(m$mrk2, levels=c("ssu", "its", "lsu"))
        ggplot(subset(m, size>1 & !is.na(value) & variable!="domain")) + geom_bar(aes(x=variable, fill=value)) + facet_grid(mrk1~mrk2) + scale_fill_manual(values=c("grey30", "red", "blue")) + labs(y="OTU count", x="taxonomic level") + theme_bw()

        #ggplot(subset(m, size>1 & !is.na(value) & variable!="domain")) + geom_bar(aes(x=variable, fill=value)) + facet_grid(mrk1~mrk2, switch="y") + scale_fill_manual(values=c("grey30", "red", "blue")) + labs(y="OTU count", x="taxonomic level") + scale_y_continuous(position = "right") + theme_bw()

        ggsave("{output.bars}", width=16, height=10)


        #proportion of comparable (assigned in both markers) OTUs
        mrk = levels(m$mrk1)
        tax = levels(m$variable)
        
        a=aggregate(oId~mrk1+mrk2+variable+value, m[m$size>1,], length)
        a$prop = NA
        a$pAss = NA
        for (m1 in mrk) {{
            for (m2 in setdiff(mrk, c(m1))) {{
                for (v1 in tax) {{
                    total = sum(subset(a, mrk1==m1 & mrk2==m2 & variable==v1)$oId)
                    ass = sum(subset(a, mrk1==m1 & mrk2==m2 & variable==v1 & (value=="same" | value=="diff"))$oId)
                    for (v2 in unique(subset(a, mrk1==m1 & mrk2==m2 & a$variable==v1)$value)) {{
                        a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$prop = a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$oId/total
                    }}
                    for (v2 in setdiff(unique(subset(a, mrk1==m1 & mrk2==m2 & a$variable==v1)$value), c("blank"))) {{
                        a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$pAss = a[a$mrk1==m1 & a$mrk2==m2 & a$variable==v1 & a$value==v2,]$oId/ass
                    }}
                }}
            }}
        }}

        #ggplot(a) + geom_line(aes(x=variable, y=pAss, group=value, color=value)) + facet_grid(mrk1~mrk2, switch="y") + scale_color_manual(values=c("grey30", "red", "blue")) + labs(y="proportion of comparable OTUs", x="taxonomic level") + theme_bw()
        a$comp=paste(a$mrk1, a$mrk2, sep="-")
        ggplot(subset(a, value=="same" & (comp %in% c("ssu-its", "ssu-lsu", "its-lsu")) & variable!="domain")) + geom_point(aes(x=variable, y=pAss, color=comp)) + geom_line(aes(x=variable, y=pAss, group=comp, color=comp), linetype=2) + theme_bw() + labs(y="proportion of comparable OTUs that have the same classification", x="taxonomic level")
        ggsave("{output.prop}", width=16, height=10)
""")

def plotChimeraEnvInput(wildcards):
    """determine input data for plotChimera rule according to sample wildcard"""
    if wildcards.sample == "all":
        return ["denovoChimera/%s.chimeraReport.tsv" % s for s in samples]
    elif wildcards.sample in isolates:
        return ["denovoChimera/%s.chimeraReport.tsv" % s for s in isoaltes[wildcards.sample]]
    else:
        return "denovoChimera/%s.chimeraReport.tsv" % wildcards.sample

rule plotChimeraEnv:
    """Plot chimera data for environmental samples"""
    input: plotChimeraEnvInput
    output: tab="denovoChimera/{sample}_chimeraTable.tsv", relative="denovoChimera/{sample}_chimerasRelativeBarplot.svg"
    run:
        R("""
        library(ggplot2)
        d=numeric(0)
        for (inFile in strsplit("{input}", " ", fixed=T)[[1]]) {{
            sample=strsplit(strsplit(inFile, "/")[[1]][2], ".", fixed=T)[[1]][1]
            n=read.table(inFile, as.is=c(2))
            colnames(n) = c("score", "query", "parentA", "parentB", "topParent", "idQM", "idQA", "idQB", "idAB", "idQT", "LY", "LN", "LA", "RY", "RN","RA", "div", "YN")
            n$YN=factor(n$YN, levels=c("Y", "?", "N"))
            n$sample=sample
            n$size = as.numeric(sub(";", "", sapply(strsplit(n$query, "="), function(x) x[4])))
            d = rbind(d, n)
        }}

        cbPalette <- c("red", "darkgrey", "blue")

        totals = aggregate(size~sample, d, sum)

        a=aggregate(size~sample+YN, d, sum)
        colnames(a) = c("sample", "YN", "count")
        a$proportion = NA
        for (i in 1:length(a$count)) {{
            a$proportion[i] = a$count[i] / totals[totals$sample==a$sample[i],]$size
        }}

        write.table(a, "{output.tab}", sep="\t", row.names=F)

        ggplot(a) + geom_bar(aes(sample, proportion, fill=YN), stat="identity") + scale_y_continuous(labels = scales::percent) + labs(x="", y="proportion of reads", title="Chimera classification in different samples") + scale_fill_manual(values=cbPalette, name=NULL, breaks=c("Y", "N", "?"), labels=c("Chimeric", "Non-chimeric", "Unclear")) + theme_bw()
        ggsave("{output.relative}", width=16, height=10)

        """)

rule combineIsolateCls:
    """De novo classify isolates for evaluation"""
    input: expand("taxonomy/{spec}_97_comb.class.tsv", spec=isolates.keys())
    output: "taxonomy/isolates_97_comb.class.tsv"
    run:
        with open(output[0], "w") as out:
            for inputFile in input:
                spec = inputFile.split("/", 1)[1].split("_", 1)[0]
                for line in open(inputFile):
                    out.write("%s\t%s" % (spec, line))
            

rule isolatePreclusterSize:
    """Generate table of pre-cluster sizes"""
    input: expand("consensus/{sample}_consensus.fasta", sample=[item for sublist in isolates.values() for item in sublist])
    output: "isoPreclusterSize.tsv"
    run:
        samp2iso = {}
        for iso, sampleList in isolates.items():
            for s in sampleList:
                samp2iso[s] = iso
        with open(output[0], "w") as out:
            for inFilePath in input:
                sample = inFilePath.split("/")[1].split("_")[0]
                iso = samp2iso[sample]
                for rec in SeqIO.parse(open(inFilePath), "fasta"):
                    size = rec.id.strip(";").rsplit(";", 1)[1].split("=")[1]
                    out.write("%s\t%s\t%s\t%s\n" % (iso, sample, rec.id, size))

################### helper functions #############################

def findTiling(hsp):
    """Find a maximum scoring subset of local hits
    
    This is done by generating a directed graph. Every match is represented by a
    Start (S...) and End (E...) node which are connected by an edge with the 
    score as weight in both directions. Compatible (in the right direction, and 
    non-overlapping) edges are connected with zero weight edges. A maximum score
    path trough the graph than represents the set of maximum scoring 
    non-overlapping matches.
    """
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
    while pre["START"][current] != start:
        path.append(pre["START"][current])
        current = pre["START"][current]
    
    return [int(p[1:]) for p in path[::-2]]

def lca(lineageStrings, stringency=1.0, 
        unidentified=["unidentified", "unclassified", "unknown"],
        ignoreIncertaeSedis=True, sizes=None):
    """Find lowest common ancestor
    
    Input is a set of strings representing classifications as path from the root
    node to the most specific available classification. Stringency define what 
    proportion of classifications need to concur at a certain level to be 
    accepted. unidetnfied gives names that should be ignored. ignoreIncertaSedis
    flag can be set to not accept incerta sedis classification. sizes can give
    wights to each classification.
    """
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



