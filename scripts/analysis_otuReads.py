import pickle

preClusterReads={}
for preInfo in snakemake.input[1:]:
    for line in open(preInfo):
        clu, seq = line.strip().split("\t")
        try:
            preClusterReads[clu].append(seq)
        except KeyError:
            preClusterReads[clu] = [seq]

clusterReads = {}
otuPreClusters = {}
repSeq = {}
for line in open(snakemake.input[0]):
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
with open(snakemake.output.size, "w") as sOut:
    for otu, readList in clusterReads.items():
        for read in readList:
            otuInf[read] = otu
        sOut.write("%s\t%i\n" % (otu, len(readList)))
with open(snakemake.output.info, "wb") as pOut:
    pickle.dump(otuInf, pOut)
with open(snakemake.output.preInfo, "wb") as eOut:
    pickle.dump(otuPreClusters, eOut)
with open(snakemake.output.repInfo, "wb") as rOut:
    pickle.dump(repSeq, rOut)
