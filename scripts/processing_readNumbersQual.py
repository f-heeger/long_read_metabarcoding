import gzip

from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    #quality filter
    for qualFileName in snakemake.input.qual:
        i=0
        with gzip.open(qualFileName, "rt") as qualFile:
            iter = SeqIO.parse(qualFile, "fastq")
            while True:
                try:
                    next(iter)
                except StopIteration:
                    break
                i+=1
        sample = qualFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
        out.write("qualFilter\t%s\t%i\n" % (sample, i))
