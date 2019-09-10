from helpers import lca
from helpers import findTiling

maxE = snakemake.config["lsuMaxEvalue"]
topPerc = snakemake.config["lsuTopPercent"]
minIdent = snakemake.config["lsuMinIdentity"]
minCov = snakemake.config["lsuMinCoverage"]
stringency = snakemake.config["lsuLcaStringency"]

taxDict={}
for line in open(snakemake.input.tax):
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

for line in open(snakemake.input.lam, encoding="latin-1"):
    total +=1
    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.strip().split("\t")
    otuId = qseqid.split("/")[0]
    if float(evalue) > maxE:
        evalueFilter += 1
        continue
    if float(pident) < minIdent:
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
topPerc = topPerc/100.0

with open(snakemake.output[0], "w") as out, open(snakemake.log[2], "w") as tLog, open(snakemake.log[1], "w") as logTax:
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
            if totalLen/min(qLen[otuId], sLen[sId])*100 < minCov:
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
            lineage = lca(goodHits, stringency)
        out.write("%s\t%s\n" % (otuId, lineage))

with open(snakemake.log[0], "w") as logOut:
    logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
    logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, maxE))
    logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, minIdent))
    logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, minCov))
