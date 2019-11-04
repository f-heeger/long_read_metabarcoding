sDat = snakemake.config["samples"]
isoLs = [sId for sId in sDat if sDat[sId]["group"] == "isolate"]

outTab = {}
oCls = {}

with open(snakemake.input[0]) as inStream:
    header = next(inStream)
    hArr = header.strip("\n").split("\t")
    assert hArr[:5] == ["otu", "totalSize", "ssuTax", "itsTax", "lsuTax"]
    samples = hArr[5:]
    iSamplesIdx = [i for i in range(len(hArr)-5) if samples[i] in isoLs]
    for line in inStream:
        lArr = line.strip("\n").split("\t")
        oId = lArr[0]
        lData = lArr[5:]
        oCls[oId] = lArr[2:5]
        if oId not in outTab:
            outTab[oId] = {}
        for ix in iSamplesIdx:
            if lData[ix] != "0":
                outTab[oId][samples[ix]] = int(lData[ix])

with open(snakemake.output[0], "wt") as out:
    out.write("otuId\tssuTax\titsTax\tlsuTax\tsampleId\treadCount\n")
    for oId, oDat in outTab.items():
        for sId, count in oDat.items():
            out.write("%s\t%s\t%s\t%i\n" % (oId, "\t".join(oCls[oId]), sId, count))

