glomeroOtu = {}
oCount = {}

with open(snakemake.input[0]) as inStream:
    headerStr = next(inStream)
    header = headerStr.strip("\n").split("\t")
    sampleIds = header[5:]
    for line in inStream:
        oId, size, ssu, its, lsu, countStr = line.strip("\n").split("\t", 5)
        if int(size) < 2:
            continue
        ssuLin = ssu.split(";")
        itsLin = its.split(";")
        lsuLin = lsu.split(";")
        if len(ssuLin) >= 4 and ssuLin[3] == "Glomeromycetes":
            try:
                glomeroOtu[oId][0] = True
            except KeyError:
                glomeroOtu[oId] = [True, False, False]
        if len(itsLin) >= 2 and itsLin[1] == "p__Glomeromycota":
            try:
                glomeroOtu[oId][1] = True
            except KeyError:
                glomeroOtu[oId] = [False, True, False]
        if len(lsuLin) >= 3 and lsuLin[2] == "Glomeromycota":
            try:
                glomeroOtu[oId][2] = True
            except KeyError:
                glomeroOtu[oId] = [False, False, True]
        oCount[oId] = dict(zip(sampleIds, [int(c) for c in countStr.split("\t")]))

with open(snakemake.output[0], "w") as out:
    out.write("sample\tglomeroReads\ttotalReads\tglomeroOtus\n")
    for sId in snakemake.config["samples"]:
        if sId not in sampleIds:
            continue
        total = 0
        glomeroLs = []
        for oId, count in oCount.items():
            if count[sId] > 0:
                total += count[sId]
                if oId in glomeroOtu:
                    glomeroLs.append(count[sId])
        out.write("%s\t%i\t%i\t%i\n" % (sId, sum(glomeroLs), total, len(glomeroLs)))
