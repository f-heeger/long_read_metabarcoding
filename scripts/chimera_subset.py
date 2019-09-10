import random

from Bio import SeqIO

random.seed("23/42")

chimera = {}
for line in open(snakemake.input.tsv):
    arr=line.strip().rsplit("\t")
    chimera[arr[1]] = arr[-1]
reads = list(SeqIO.parse(open(snakemake.input.seqs), "fasta"))
subset = random.sample(reads, snakemake.params.N)

with open(snakemake.output[0], "w") as out:
    for r, read in enumerate(subset):
        read.description = chimera[read.id]
        read.id = "%i_%s" % (r, read.id)
        out.write(read.format("fasta"))

