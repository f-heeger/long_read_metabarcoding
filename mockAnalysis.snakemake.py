import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

import networkx as nx

from snakemake.utils import min_version, R

shell.prefix("sleep 10; ") #work around to desl with "too quck" rule execution and slow samba

configfile: "config.json"

rule all:
    input: "mock/clusterGraph/Lib4-0018_clusterGraphCls.tsv", "mock/clusterGraph/Lib4-0018_clusterGraphEdges.tsv", expand("mock/clusterAln/Lib4-0018_compAln_0_{marker}.fasta", marker=["SSU", "ITS1", "58S", "ITS2", "LSU"]),

rule removeChimera:
    input: seqs="primers/{sample}_minLen.fasta", ref="isolateSeqs.fasta"
    output: fasta="mock/refChimera/{sample}.nochimera.fasta", tsv="mock/refChimera/{sample}.chimeraReport.tsv"
    log: "mock/logs/{sample}_refChimera.log"
    threads: 6
    shell:
        "%(vsearch)s --uchime_ref {input.seqs} --db {input.ref} --nonchimeras {output.fasta} --uchimeout {output.tsv} --threads {threads} &> {log}" % config

rule itsx:
    input: "mock/refChimera/{sample}.nochimera.fasta"
    output: "mock/itsx/{sample}.SSU.fasta", "mock/itsx/{sample}.ITS1.fasta", "mock/itsx/{sample}.5_8S.fasta", "mock/itsx/{sample}.ITS2.fasta", "mock/itsx/{sample}.LSU.fasta", "mock/itsx/{sample}.summary.txt", "mock/itsx/{sample}.positions.txt", "mock/itsx/{sample}.full.fasta"
    threads: 6
    log: "mock/logs/{sample}_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o mock/itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 500 -E 1e-4 2> {log}" % config
        

rule indiDerep:
    input: "mock/itsx/{sample}.{marker}.fasta"
    output: fasta="mock/indiDerep/{sample}_{marker}.derep.fasta", info="mock/indiDerep/{sample}_{marker}.repseq.pic", txt="mock/indiDerep/{sample}_{marker}.uc.txt"
    log: "mock/logs/{sample}_{marker}_derep.log"
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
    input: "mock/indiDerep/{sample}_{marker}.derep.fasta"
    output: fasta="mock/indiCluster/{sample}_{marker}_clu.fasta", uc="mock/indiCluster/{sample}_{marker}_clu.uc.tsv"
    log: "mock/logs/{sample}_indiClustering_{marker}.log"
    threads: 6
    shell:
        "%(vsearch)s --cluster_size {input} --iddef 0 --id 0.97 --minsl 0.9 --sizein --sizeout --centroids {output.fasta} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

rule cp58s:
    input: "mock/itsx/{sample}.5_8S.fasta"
    output: "mock/itsx/{sample}.58S.fasta"
    shell:
        "cp {input} {output}"

rule indiCluReads:
    input: clsInfo="mock/indiCluster/{sample}_{marker}_clu.uc.tsv", repseq="mock/indiDerep/{sample}_{marker}.repseq.pic"
    output: info="mock/indiCluster/{sample}_{marker}_cluInfo.pic"
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
    input: ssu="mock/indiCluster/{sample}_SSU_cluInfo.pic", its1="mock/indiCluster/{sample}_ITS1_cluInfo.pic", r58s="mock/indiCluster/{sample}_58S_cluInfo.pic", its2="mock/indiCluster/{sample}_ITS2_cluInfo.pic", lsu="mock/indiCluster/{sample}_LSU_cluInfo.pic"
    output: nodes="mock/clusterGraph/{sample}_clusterGraphNodes.tsv", edges="mock/clusterGraph/{sample}_clusterGraphEdges.tsv"
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
    input: edges="mock/clusterGraph/{sample}_clusterGraphEdges.tsv"
    output: comp="mock/clusterGraph/{sample}_components.txt"
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
    input: ssu="mock/indiCluster/{sample}_SSU_cluInfo.pic", its1="mock/indiCluster/{sample}_ITS1_cluInfo.pic", r58s="mock/indiCluster/{sample}_58S_cluInfo.pic", its2="mock/indiCluster/{sample}_ITS2_cluInfo.pic", lsu="mock/indiCluster/{sample}_LSU_cluInfo.pic", cls="mapping/assignment/{sample}_assignments.tsv", comp="mock/clusterGraph/{sample}_components.txt"
    output: tab="mock/clusterGraph/{sample}_components_class.tsv", lab="mock/clusterGraph/{sample}_clusterGraphCls.tsv"
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
    input: comp="mock/clusterGraph/{sample}_components.txt", fasta="mock/itsx/{sample}.{marker}.fasta", info="mock/indiCluster/{sample}_{marker}_cluInfo.pic"
    output: dynamic("mock/clusterAln/{sample}_compAln_{comp}_{marker}.fasta")
    log: comp="mock/logs/{sample}_components.log", aln="mock/logs/{sample}_cluAlign.log"
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
                
                baseFileName = "mock/clusterAln/tmp_comp%i_%s" % (l, wildcards.marker)
                #write all reads in the component to temp file
                readFilePath = "%s_reads.fasta" % baseFileName
                alignmentPath="mock/clusterAln/%s_compAln_%i_%s.fasta" % (wildcards.sample, l, wildcards.marker)
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
                shell("rm mock/clusterAln/tmp_*")

