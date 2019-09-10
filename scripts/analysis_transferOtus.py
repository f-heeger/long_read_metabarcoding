import pickle

from Bio import SeqIO

repSeq = pickle.load(open(snakemake.input[0], "rb"))
rep2otu = dict(zip(repSeq.values(), repSeq.keys()))
with open(snakemake.output[0], "w") as out:
    for rec in SeqIO.parse(open(snakemake.input[1]), "fasta"):
        readId = rec.id.split("=", 1)[1].split(";",1)[0]
        if readId in rep2otu:
            rec.id = "%s/%s" % (rep2otu[readId], snakemake.wildcards.marker)
            out.write(rec.format("fasta"))
