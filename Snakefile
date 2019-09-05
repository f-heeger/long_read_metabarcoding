import pickle
import json
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from snakemake.utils import min_version, R



shell.prefix("sleep 10; ") #work around to deal with "too quick" rule execution and slow samba

configfile: "config.json"

sampleInfo = json.load(open("samples.json"))

samples = [sId for sId, sType in sampleInfo["sampleType"].items() if sType == "env"]
isolates = {}
for sId, sType in sampleInfo["sampleType"].items():
    if sType == "isolate":
        try:
            isolates[sampleInfo["sampleName"][sId]].append(sId)
        except KeyError:
            isolates[sampleInfo["sampleName"][sId]] = [sId]

allSamples = list(sampleInfo["sampleName"].keys())

sampleName = sampleInfo["sampleName"]

mockSamples = [sId for sId, sType in sampleInfo["sampleType"].items() if sType == "mock"]

####################################################################
# includes

include: "rules/createDBs.snakemake.py"
include: "rules/readProcessing.snakemake.py"
include: "rules/mapping.snakemake.py"
include: "rules/chimera_analysis.snakemake.py"
include: "rules/mockAnalysis.snakemake.py"
include: "rules/analysis.snakemake.py"

rule all:
    input: "all_otu97_table.tsv", "taxonomy/all_97_comb.class.tsv", "all_clsComp_depth.svg", "all_clsComp_depth_fungi.svg", "all_clsComp_basic.svg", "all_clsDiffStat.svg", "taxonomy/Lib4-0018_97_combToCorr.class.tsv", "chimeraCyclesRelativeBarplot.svg", "chimera_comp_sankey.svg", expand(["mapping/{stage}MockComp.svg", "mapping/{stage}ErrorRates.svg"], stage=["raw", "filtered"]), "readNumbers.svg", "mock/clusterGraph/Lib4-0018_clusterGraphCls.tsv", "mock/clusterGraph/Lib4-0018_clusterGraphEdges.tsv", "mock/clusterGraph/Lib4-0018_clusterGraphClsLab.tsv"

rule metabarcoding:
    input: "all_otu97_table.tsv", "all_clsComp_depth.svg"
