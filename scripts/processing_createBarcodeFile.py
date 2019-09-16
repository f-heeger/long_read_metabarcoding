bc = {}
for sId, sample in snakemake.config["samples"].items():
    bc[sample["fwdBarcodeId"]] = sample["fwdBarcodeSeq"]
    bc[sample["revBarcodeId"]] = sample["revBarcodeSeq"]

with open(snakemake.output[0], "wt") as out:
    for bcIt in bc.items():
        out.write(">%s\n%s\n" % bcIt)
