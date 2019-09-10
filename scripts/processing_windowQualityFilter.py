import gzip

from Bio import SeqIO

from helpers import homopoly

winSize = snakemake.config["qualityWindowSize"]
minQual = snakemake.config["windowMinQuality"]
removed = 0
with open(snakemake.output.good, "wt") as out, open(snakemake.output.stat, "w") as statOut:
    for read in SeqIO.parse(gzip.open(snakemake.input.fastq, "rt"), "fastq"):
        anyRemoved = False
        for i in range(len(read)-winSize):
            tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"][i:i+winSize] ]) / winSize
            tQual = 1 - tError
            tRemoved=False
            if (tQual) < minQual:
                tRemoved = True
                anyRemoved = True
            kmer=str(read.seq[i:i+winSize])
            hp = homopoly(kmer)
            statOut.write("%s\t%s\t%i\t%i\t%f\t%s\t%s\t%i\t%i\t%i\t%i\t%i\n" % (snakemake.wildcards.sample, read.id, i, len(read), tQual, tRemoved, kmer, kmer.count("A"), kmer.count("C"), kmer.count("G"), kmer.count("T"), hp))
        if not anyRemoved:
            out.write(read.format("fasta"))
        else:
            removed += 1

with open(snakemake.log[0], "w") as logFile:
    logFile.write("%s: %i reads removed because in a window of size %i quality droped below %f\n" % (snakemake.wildcards.sample, removed, winSize, minQual))
