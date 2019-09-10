from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    for read in SeqIO.parse(open(snakemake.input[0]), "fastq"):
        out.write(read.format("fasta"))
