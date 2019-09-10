from Bio import SeqIO


file2sample = dict(zip([s["path"] for s in config["samples"].values()],
                       config["samples"].keys()))

with open(snakemake.output[0], "w") as out:
    #raw
    for rawFileName in snakemake.input.raw:
        i=0
        print(rawFileName)
        import pdb; pdb.set_trace()
        with open(rawFileName) as rawFile:
            iter = SeqIO.parse(rawFile, "fastq")
            while True:
                try:
                    next(iter)
                except StopIteration:
                    break
                i+=1
            
            out.write("raw\t%s\t%i\n" % (file2sample[rawFileName], i))
