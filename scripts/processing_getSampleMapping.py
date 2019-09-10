import pickle

from Bio import SeqIO

preClu2sample = {}
for inputFile in snakemake.input:
    sample = inputFile.rsplit("/", 1)[-1].split(".", 1)[0]
    for rec in SeqIO.parse(open(inputFile), "fasta"):
        preClu2sample[rec.id] = sample
pickle.dump(preClu2sample, open(snakemake.output.sample, "wb"))
