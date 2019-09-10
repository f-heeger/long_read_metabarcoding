from Bio import SeqIO

seqQ = {}
pinQ = {}
with open(snakemake.input.qual) as qual:
    header = next(qual)
    for line in qual:
        arr = line.strip().split("\t")
        seqQ[arr[0]] = float(arr[5])
        try:
            pinQ[arr[0]] = float(arr[12])
        except ValueError:
            if arr[12] == "NULL":
                pinQ[arr[0]] = None
            else:
                raise

with open(snakemake.output[0], "w") as out:
    total=0
    badS=0
    badP=0
    for rec in SeqIO.parse(snakemake.open(input.fasta), "fasta"):
        total += 1
        acc = rec.id.split(".")[0]
        if seqQ[acc] < snakemake.params.minSeqQ:
            badS += 1
            continue
        if pinQ[acc] is None or pinQ[acc] < snakemake.params.minPinQ:
            badP += 1
            continue
        out.write(rec.format("fasta"))
    open(snakemake.log[0], "w").write("%i sequences read. %i sequences removed because of sequence quality < %f.\n%i sequences removed because of pintail quality < %f.\n" % (total, badS, snakemake.params.minSeqQ, badP, snakemake.params.minPinQ))
