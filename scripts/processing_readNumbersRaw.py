import gzip

from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    #raw
    for rawFileName in snakemake.input:
        i=0
        with gzip.open(rawFileName, "rt") as rawFile:
            iter = SeqIO.parse(rawFile, "fastq")
            while True:
                try:
                    next(iter)
                except StopIteration:
                    break
                i+=1
            sample = rawFileName.split("/")[-1].split(".")[0]
            out.write("raw\t%s\t%i\n" % (sample, i))
