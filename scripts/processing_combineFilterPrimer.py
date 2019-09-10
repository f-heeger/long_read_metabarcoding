from Bio import SeqIO

with open(snakemake.output.fastq, "wt") as out:
    for rec in SeqIO.parse(snakemake.input.fwd,"fasta"):
        out.write(rec.format("fasta"))
    for rec in SeqIO.parse(snakemake.input.rev,"fasta"):
        newRec = rec.reverse_complement()
        newRec.id = "%s/rc" % rec.id
        newRec.description = ""
        out.write(newRec.format("fasta"))

