import gzip

from Bio import SeqIO

maxLen = snakemake.config["maxReadLen"]
minLen = snakemake.config["minReadLen"]
nrLong = 0
nrShort = 0
with gzip.open(snakemake.output.right, "wt") as out, gzip.open(snakemake.output.long, "wt") as tooLong, gzip.open(snakemake.output.short, "wt") as tooShort:
    for read in SeqIO.parse(gzip.open(snakemake.input[0], "rt"), "fastq"):
        if len(read) <= minLen:
            tooShort.write(read.format("fastq"))
            nrShort += 1
        elif len(read) <= maxLen:
            out.write(read.format("fastq"))
        else:
            tooLong.write(read.format("fastq"))
            nrLong += 1
with open(snakemake.log[0], "w") as logFile:
    logFile.write("%i reads were removed because they were longer than %i\n" % (nrLong, maxLen))
    logFile.write("%i reads were removed because they were shorter than %i\n" % (nrShort, minLen))
