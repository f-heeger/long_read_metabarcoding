import gzip

from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    #max lenght filter
    for lenFileName in snakemake.input.length:
        i=0
        with gzip.open(lenFileName, "rt") as lenFile:
            iter = SeqIO.parse(lenFile, "fastq")
            while True:
                try:
                    next(iter)
                except StopIteration:
                    break
                i+=1
        sample = lenFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
        out.write("lenFilter\t%s\t%i\n" % (sample, i))
