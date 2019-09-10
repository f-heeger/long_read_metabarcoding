from Bio import SeqIO

accRank = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
rank = {}
for line in open(snakemake.input.tax):
    path, tId, tRank, remark, release = line.strip("\n").split("\t")
    rank[path.strip(";").replace(" ", "_")] = tRank
with open(snakemake.output.tax, "w") as tOut:
    for rec in SeqIO.parse(open(snakemake.input.fasta), "fasta"):
        tax = rec.description.split(" ", 1)[1].replace(" ", "_")
        taxArr = tax.strip(";").split(";")
        newTax = dict(zip(accRank, ["Incertae sedis"]*8))
        for i in range(len(taxArr)-1):
            try:
                tRank = rank[";".join(taxArr[:i+1])]
            except KeyError:
                if ";".join(taxArr[:i+1]) == "Bacteria;RsaHf231":
                    tRank="phylum" #work around for error in input file (v128), also mentioned here: http://blog.mothur.org/2017/03/22/SILVA-v128-reference-files/
                else:
                    raise
            if tRank in accRank:
                newTax[tRank] = taxArr[i]
        newTax["species"] = taxArr[-1]
        tOut.write("%s\t%s;\n" % (rec.id, ";".join([newTax[r] for r in accRank])))
