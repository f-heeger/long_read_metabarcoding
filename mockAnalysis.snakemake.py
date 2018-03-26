

rule indiItsx:
    """Run ITSx on each non-chimeric pre-cluster """
    input: "refChimera/{sample}.nochimera.fasta"
    output: "mock/itsx/{sample}.SSU.fasta", "mock/itsx/{sample}.ITS1.fasta", "mock/itsx/{sample}.5_8S.fasta", "mock/itsx/{sample}.ITS2.fasta", "mock/itsx/{sample}.LSU.fasta", "mock/itsx/{sample}.summary.txt", "mock/itsx/{sample}.positions.txt", "mock/itsx/{sample}.full.fasta"
    threads: 6
    log: "mock/logs/{sample}_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o mock/itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 500 -E 1e-4 2> {log}" % config
        

rule indiDerep:
    """Dereplicate each marker independently"""
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
    """Cluster each marker independently"""
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
    """Collect information which read is in which cluster"""
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
    """Check which reads are in one cluster and which reads are shared between clusters.
    
       Generate information for edges (shared reads between clusters of diferent markers) and nodes (read per cluster) in the cluster graph.
    """
    input: ssu="mock/indiCluster/{sample}_SSU_cluInfo.pic", its1="mock/indiCluster/{sample}_ITS1_cluInfo.pic", r58s="mock/indiCluster/{sample}_58S_cluInfo.pic", its2="mock/indiCluster/{sample}_ITS2_cluInfo.pic", lsu="mock/indiCluster/{sample}_LSU_cluInfo.pic"
    output: nodes="mock/clusterGraph/{sample}_clusterGraphNodes.tsv", edges="mock/clusterGraph/{sample}_clusterGraphEdges.tsv"
    params: minCluSize=3, minCluOverlap=3
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
            out.write("cluster\ttype\tsize\tlogSize\n")
            for cId, reads in nodes.items():
                out.write("%s\t%s\t%i\t%f\n" % (cId, cId.split("|")[2], len(reads), math.log(len(reads), 2)))
                    
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
    """Find connected components in the cluster graph"""
    input: edges="mock/clusterGraph/{sample}_clusterGraphEdges.tsv", nodes="mock/clusterGraph/{sample}_clusterGraphNodes.tsv"
    output: comp="mock/clusterGraph/{sample}_components.txt"
    run:
        G=nx.Graph()
        with open(input.nodes) as nodeStream:
            _ = next(nodeStream) # header
            for line in nodeStream:
                nId, _ = line.strip().split("\t", 1)
                G.add_node(nId)
        with open(input.edges) as edgeStream:
            _ = next(edgeStream) # header
            for line in edgeStream:
                start, end, weight = line.strip().split("\t")
                G.add_edge(start, end)
        with open(output.comp, "w") as out:
            for comp in nx.connected_components(G):
                out.write("\t".join(comp)+"\n")

rule indiCluClass:
    """Find classifications for each read to label the cluster graph with"""
    input: ssu="mock/indiCluster/{sample}_SSU_cluInfo.pic", its1="mock/indiCluster/{sample}_ITS1_cluInfo.pic", r58s="mock/indiCluster/{sample}_58S_cluInfo.pic", its2="mock/indiCluster/{sample}_ITS2_cluInfo.pic", lsu="mock/indiCluster/{sample}_LSU_cluInfo.pic", cls="mapping/assignment/{sample}_filtered_assignments.tsv", comp="mock/clusterGraph/{sample}_components.txt"
    output: tab="mock/clusterGraph/{sample}_components_class.tsv", lab="mock/clusterGraph/{sample}_clusterGraphClsLab.tsv", cls="mock/clusterGraph/{sample}_clusterGraphCls.tsv"
    run:
        reads = {
        "SSU": pickle.load(open(input.ssu, "rb")),
        "ITS1": pickle.load(open(input.its1, "rb")),
        "5.8S": pickle.load(open(input.r58s, "rb")),
        "ITS2": pickle.load(open(input.its2, "rb")),
        "LSU": pickle.load(open(input.lsu, "rb"))
        }
        readCls = {}
        uCls = set([])
        for line in open(input.cls):
            read, cls = line.strip().split("\t")[:2]
            readCls[read] = cls
            uCls.add(cls)
        clsOrder = list(uCls)
        with open(output.tab, "w") as out, open(output.lab, "w") as labOut, open(output.cls, "w") as clsOut:
            labOut.write("nodeName\tclassification\n")
            clsOut.write("nodeName\t%s\n" % "\t".join(clsOrder))
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
                    clsOut.write("%s\t%s\n" % (clu, "\t".join([str(cluCls.get(c, 0)) for c in clsOrder])))
