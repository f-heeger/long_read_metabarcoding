from herlpers import lca

maxE = snakemake.config["itsMaxEvalue"]
topPerc = snakemake.config["itsTopPercent"]
minIdent = snakemake.config["itsMinIdentity"]
minCov = snakemake.config["itsMinCoverage"]
stringency = snakemake.config["itsLcaStringency"]

taxDict = {}
for line in open(snakemake.input.tax):
    sh, tax = line.strip().split("\t")
    taxDict[sh] = tax
logOut = open(snakemake.log[0], "w")
logTax = open(snakemake.log[1], "w")
classifi = {}
seqLength = {}
seqNr = 0
total = 0
evalueFilter = 0
identFilter = 0
covFilter = 0
for line in open(snakemake.input.lam, encoding="latin-1"):
    total +=1
    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.strip().split("\t")
    readId = qseqid.split("|")[0]
    if float(evalue) > maxE:
        evalueFilter += 1
        continue
    if float(pident) < minIdent:
        identFilter +=1
        continue
    mLen = min(int(qlen), int(slen))
    if float(length)/mLen*100 < minCov:
        covFilter += 1
        continue
    linStr = taxDict[sseqid]
    try:
        classifi[readId].append((linStr, float(bitscore)))
    except KeyError:
        classifi[readId] = [(linStr, float(bitscore))]
logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, maxE))
logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, minIdent))
logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, minCov))
topPerc = topPerc/100.0
with open(snakemake.output[0], "w") as out:
    for key, hits in classifi.items():
        sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
        cutoff = 0
        while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
            cutoff += 1
        goodHits = [hit[0] for hit in sortedHits[:cutoff]]
        for h in goodHits:
            logTax.write("%s\t%s\n" % (key, h))
        lineage = lca(goodHits, stringency)
        out.write("%s\t%s\n" % (key, lineage))
try:
    logOut.close()
except:
    pass
try:
    logTax.close()
except:
    pass
