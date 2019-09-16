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
    sId, sName, fwdBcId, fwdBcSeq, revBcId, revBcSeq = line.strip("\n").split("\t")
    samples[sId] = {"name": sName, 
                    "fwdBarcodeId": fwdBcId, "revBarcodeId": revBcId,
                    "fwdBarcodeSeq": fwdBcSeq, "revBarcodeSeq": revBcSeq,
                    }

config["samples"] = samples

####################################################################
# includes

include: "rules/createDBs.snakemake.py"
include: "rules/readProcessing.snakemake.py"
#include: "rules/mapping.snakemake.py"
include: "rules/chimera_analysis.snakemake.py"
#include: "rules/mockAnalysis.snakemake.py"
include: "rules/analysis.snakemake.py"

rule all:
    input: "all_otu97_table.tsv", "taxonomy/all_97_comb.class.tsv", "all_clsComp_depth.svg", "all_clsComp_depth_fungi.svg", "all_clsComp_basic.svg", "all_clsDiffStat.svg", "taxonomy/Lib4-0018_97_combToCorr.class.tsv", "chimeraCyclesRelativeBarplot.svg", "chimera_comp_sankey.svg", expand(["mapping/{stage}MockComp.svg", "mapping/{stage}ErrorRates.svg"], stage=["raw", "filtered"]), "readNumbers.svg", "mock/clusterGraph/Lib4-0018_clusterGraphCls.tsv", "mock/clusterGraph/Lib4-0018_clusterGraphEdges.tsv", "mock/clusterGraph/Lib4-0018_clusterGraphClsLab.tsv"

rule metabarcoding:
    input: "all_otu97_table.tsv", "all_clsComp_depth.pdf"
