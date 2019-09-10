from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    #primer filter
    for primerFileName in snakemake.input.primer:
        i=0
        with open(primerFileName) as primerFile:
            iter = SeqIO.parse(primerFile, "fasta")
            while True:
                try:
                    next(iter)
                except StopIteration:
                    break
                i+=1
        sample = primerFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
        out.write("primerFilter\t%s\t%i\n" % (sample, i))
