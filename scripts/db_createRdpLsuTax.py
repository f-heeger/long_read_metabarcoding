from Bio import SeqIO

accRank = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
with open(snakemake.output[0], "w") as out:
    for rec in SeqIO.parse(open(snakemake.input[0]), "fasta"):
        seqIdStr, linStr = rec.description.strip().split("\t")
        linArr = linStr.split("=", 1)[1].split(";")
        spec = seqIdStr.split(";", 1)[0].split(" ", 1)[1]
        newTax = dict(zip(accRank, [None]*8))
        for i in range(0, len(linArr), 2):
            if linArr[i+1] in accRank:
                newTax[linArr[i+1]] = linArr[i]
        newTax["kingdom"] = "Fungi"
        newTax["domain"] = "Eukaryota"
        newTax["species"] = spec
        #########################################
        #work arounds for wired cases
#                if spec == "Rhizophlyctis rosea":
#                    newTax = {"kingdom": "Fungi", "domain": "Eukaryota", "pyhlum": "Chytridiomycota", "class": "Chytridiomycetes", "order": "Spizellomycetales", "family": "Spizellomycetaceae", "genus": "Rhizophlyctis"}
        #########################################
        found=False
        for rank in accRank[-2::-1]:
            if newTax[rank] is None:
                if found:
                    newTax[rank] = "Incertae sedis"
                else:
                    newTax[rank] = "unclassified"
            
        out.write("%s\t%s;\n" % (rec.id, ";".join([newTax[r] for r in accRank])))

