from Bio import SeqIO

file2sample = dict(zip([s["path"] for s in snakemake.config["samples"].values()],
                       snakemake.config["samples"].keys()))

with open(snakemake.output[0], "w") as out:
    #raw
    for rawFileName in snakemake.input:
        i=0
        with open(rawFileName) as rawFile:
            iter = SeqIO.parse(rawFile, "fastq")
            while True:
                try:
                    next(iter)
                except StopIteration:
                    break
                i+=1
            
            out.write("raw\t%s\t%i\n" % (file2sample[rawFileName], i))
