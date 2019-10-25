from snakemake.utils import min_version, R



shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow samba

configfile: "config.json"

comp = {"A": "T", "T": "A", "C": "G", "G": "C",
        "W": "W", "M": "K", "K": "M", "R": "Y",
        "Y": "R", "S": "S", "H": "D", "D": "H",
        "V": "B", "B": "V", "N": "N"}

config["rv_fwd_primer"] = "".join(comp[base] for base in config["fwd_primer"][::-1])
config["rv_rev_primer"] = "".join(comp[base] for base in config["rev_primer"][::-1])

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
include: "rules/chimera_analysis.snakemake.py"
include: "rules/analysis.snakemake.py"

rule all:
    input: "all_otu97_table.tsv", "all_clsComp_depth.svg", "all_clsComp_depth_fungi.svg", "all_clsDiffStat.svg", "readNumbers.svg"

