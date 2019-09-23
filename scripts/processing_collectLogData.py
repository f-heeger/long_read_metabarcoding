d = {}
#raw
for line in open(snakemake.input.raw):
    _, sample, num = line.strip("\n").split("\t")
    d[sample] = {"raw": int(num)}

#length
for inFile in snakemake.input.length:
    sample = inFile.split("/")[1].rsplit("_", 1)[0]
    with open(inFile) as inStream:
        d[sample]["tooLong"] = int(next(inStream).split(" ", 1)[0])
        d[sample]["tooShort"] = int(next(inStream).split(" ", 1)[0])
#qual
for inFile in snakemake.input.qual:
    sampleStr, number, rest = open(inFile).read().split(" ", 2)
    d[sampleStr[:-1]]["lowQual"] = int(number)
#window
for inFile in snakemake.input.window:
    sampleStr, number, rest = open(inFile).read().split(" ", 2)
    d[sampleStr[:-1]]["lowWinQual"] = int(number)
primer = {}
#primer1
for inFile in snakemake.input.primer1:
    sample = inFile.split("/")[1].rsplit("_", 2)[0]
    primer[sample] = {}
    for line in open(inFile):
        if line.startswith("Total reads processed"):
            label, data = line.split(":")
            primer[sample]["total"] = int(data.strip().replace(",",""))
        elif line.startswith("Reads with adapters"):
            label, data = line.split(":")
            primer[sample]["forward"] = int(data.strip().split(" ")[0].replace(",",""))
            break
#primer2
for inFile in snakemake.input.primer2:
    sample = inFile.split("/")[1].rsplit("_", 2)[0]
    for line in open(inFile):
        if line.startswith("Total reads processed"):
            label, data = line.split(":")
            total = int(data.strip().replace(",",""))
            assert primer[sample]["total"] == total
        elif line.startswith("Reads with adapters"):
            label, data = line.split(":")
            primer[sample]["reverse"] = int(data.strip().split(" ")[0].replace(",",""))
            break

for sId, pData in primer.items():
    d[sId]["noPrimer"] = pData["total"] - pData["forward"] - pData["reverse"]
#chimera
for inFile in snakemake.input.chimera:
    sample = inFile.split("/")[1].rsplit("_", 1)[0]
    lines = open(inFile).readlines()
    #Found 1 (0.6%) chimeras, 167 (99.4%) non-chimeras,
    chimStr, nonChimStr, _ = lines[9].split(",")
    chim = chimStr.split(" ", 2)[1]
    #and 0 (0.0%) borderline sequences in 168 unique sequences.
    susp = lines[10].split(" ")[1]
    d[sample]["chimeras"] = int(chim) + int(susp)
#itsx
for inFile in snakemake.input.itsx:
    sample = inFile.split("/")[1].rsplit(".", 2)[0]
    for line in open(inFile):
        if line.startswith("Number of sequences in input file"):
            text, totalStr = line.strip("\n").split("\t")
        if line.startswith("Sequences detected as ITS by ITSx:"):
            text, foundStr = line.strip("\n").split("\t")
            d[sample]["noITS"] = int(totalStr) - int(foundStr)

with open(snakemake.output[0], "w") as out:
    header = ["raw", "tooLong", "tooShort", "lowQual", "lowWinQual", "noPrimer", "chimeras", "noITS"]
    out.write("sample\t"+"\t".join(header) + "\n")
    for sId, data in d.items():
        out.write("%s\t%s\n" % (sId, "\t".join([str(data.get(h, "NA")) for h in header])))
    
