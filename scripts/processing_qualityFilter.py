import gzip

from Bio import SeqIO

minQual = snakemake.config["minReadQual"]
removed = 0

with gzip.open(snakemake.output.good, "wt") as out, gzip.open(snakemake.output.bad, "wt") as trash, open(snakemake.output.info, "w") as info:
    for read in SeqIO.parse(gzip.open(snakemake.input.fastq, "rt"), "fastq"):
        try:
            tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"]]) / len(read)
        except Exception:
            print(read.id)
            raise
        info.write("%s\t%f\n" % (read.id, 1-tError))
        if (1-tError) < minQual:
            removed += 1
            trash.write(read.format("fastq"))
        else:
            out.write(read.format("fastq"))
open(snakemake.log[0], "w").write("%s: %i reads removed because quality < %f\n" % (snakemake.wildcards.sample, removed, minQual))
