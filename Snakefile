from snakemake.utils import min_version, R



shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow samba

configfile: "config.json"

samples = {}
for line in open("samples.tsv"):
    if line[0] == "#":
        continue
    sId, sName, sraId, fwdBcId, fwdBcSeq, revBcId, revBcSeq, group = line.strip("\n").split("\t")
    samples[sId] = {"name": sName, "sraId": sraId,
                    "fwdBarcodeId": fwdBcId, "revBarcodeId": revBcId,
                    "fwdBarcodeSeq": fwdBcSeq, "revBarcodeSeq": revBcSeq,
                    "group": group
                    }

config["samples"] = samples

####################################################################
# includes

include: "rules/createDBs.snakemake.py"
include: "rules/readProcessing.snakemake.py"
include: "rules/analysis.snakemake.py"

rule all:
    input: "all_otu97_table.tsv", "all_clsComp_depth.svg", "all_clsComp_depth_fungi.svg", "all_clsDiffStat.svg", "readNumbers/readNumbersFiltering.svg", "readsAfterFiltering.pdf"

