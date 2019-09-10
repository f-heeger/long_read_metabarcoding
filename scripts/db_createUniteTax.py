from Bio import SeqIO

with open(snakemake.output.tax, "w") as tOut, open(snakemake.output.fasta, "w") as fOut:
    for rec in SeqIO.parse(open(snakemake.input[0]), "fasta"):
        arr = rec.id.split("|")
        sh = arr[2]
        tax = arr[-1]
        tOut.write("%s\t%s\n" % (sh, tax))
        rec.id = sh
        fOut.write(rec.format("fasta"))
