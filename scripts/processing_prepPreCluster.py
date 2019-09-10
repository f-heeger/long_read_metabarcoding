from Bio import SeqIO

qual = {}
for line in open(snakemake.input.qual):
    rId, tQual = line.strip().split("\t")
    qual[rId] = tQual

readList = []
for read in SeqIO.parse(open(snakemake.input.fasta), "fasta"):
    if read.id.endswith("rc"):
        readId = read.id.rsplit("/", 1)[0]
    else:
        readId = read.id
    readList.append((read, qual[readId]))

with open(snakemake.output[0], "w") as out:
    for read, qual in sorted(readList, key=lambda x: x[1], reverse=True):
        read.description=str(qual)
        out.write(read.format("fasta"))
