import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from snakemake.utils import min_version, R

shell.prefix("sleep 10; ") #work around to desl with "too quck" rule execution and slow samba

configfile: "config.json"

barcodes=["0009", "0075", "0018", "0095", "0027", "0034", "0056"]


samples = []
for lib, bcList in config["samples"].items():
    for bc in bcList: 
        samples.append("%s_%s" % (lib, bc))

#print(sorted(samples))

rule all:
    input: expand("taxonomy/{sample}_comb.class.tsv", sample=["Lib4_0018", "Lib3_0075", "Lib3_0034", "Lib7_0075", "Lib7_0034"])
    #"QC/multiqc_report.html", "readNumbers.pdf", expand("taxonomy/{sample}.clu.class.tsv", sample=["Lib%i_0075" % i for i in range(1,9)]+["Lib%i_0034" % i for i in [1,2,3,5,6,7,8]]), expand("clusters2/{sample}_cluster2.size.tsv", sample=samples), #"clusters/all_cluster_persample.tsv"

rule unpack:
    input: "raw/8_libs_Mar17/Ampl.Lib{cellNr}.SC1+2_barcoded-fastqs.tgz"
    output: dynamic("raw/Lib{cellNr}/{barcode}_Forward--{barcode}_Forward.fastq")
    shell:
        "mkdir -p raw/Lib{wildcards.cellNr}; tar -xzf {input} -C raw/Lib{wildcards.cellNr}; touch raw/Lib{wildcards.cellNr}/*"

rule renameRawfile:
    input: "raw/Lib{cellNr}/{barcode}_Forward--{barcode}_Forward.fastq"
    output: "raw/Lib{cellNr}_{barcode}.fastq"
    shell:
        "mv {input} {output}"

rule fastqc:
    input: "raw/{sample}.fastq"
    output: "QC/{sample}_fastqc.html"
    threads: 6
    shell: "%(fastqc)s -o QC -t {threads} {input}" % config

rule multiqc:
    input: expand("QC/{sample}_fastqc.html", sample=samples)
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    shell:
        "%(multiqc)s -f --interactive -o QC QC/*_fastqc.zip" % config

rule lengthFilter:
    input: fastq="raw/{sample}.fastq"
    output: "lenFilter/{sample}_rightLen.fastq", "lenFilter/{sample}_wrongLen.fastq"
    log: "logs/{sample}_lenFilter.log"
    params: maxLen = 6000
    run:
        toLong = 0
        with open(output[0], "w") as out, open(output[1], "w") as trash:
            for read in SeqIO.parse(open(input.fastq), "fastq"):
                if len(read) <= params.maxLen:
                    out.write(read.format("fastq"))
                else:
                    trash.write(read.format("fastq"))
                    toLong += 1
        open(log[0], "w").write("%s: %i reads were removed because they were longer than %i\n" % (wildcards.sample, toLong, params.maxLen))

rule qualityFilter:
    input: fastq="lenFilter/{sample}_rightLen.fastq"
    output: "qualFilter/{sample}_goodQual.fastq", "qualFilter/{sample}_badQual.fastq"
    params: minQual=0.99
    log: "logs/{sample}_qualityFilter.log"
    run:
        removed = 0
        with open(output[0], "w") as out, open(output[1], "w") as trash:
            for read in SeqIO.parse(open(input.fastq), "fastq"):
                rId, qual, _ = read.description.split(" ")
                if float(qual) < params.minQual:
                    removed += 1
                    trash.write(read.format("fastq"))
                else:
                    out.write(read.format("fastq"))
        open(log[0], "w").write("%s: %i reads removed because quality < %f\n" % (wildcards.sample, removed, params.minQual))

rule init_filterPrimer:
    input: "qualFilter/{sample}_goodQual.fastq"
    output: fastq="primers/{sample}_primer.fastq"
    log: "logs/{sample}_primer.log"
    threads: 3
    params: minovl=10
    run:
        shell("%(cutadapt)s -g forward=%(fwd_primer)s -g reverse=%(rev_primer)s -O {params.minovl} --trimmed-only -o primers/{wildcards.sample}_{{name}}.fastq {input} &> {log}" % config)
        config["fwd_primer_rc"] = str(Seq(config["fwd_primer"], IUPACAmbiguousDNA()).reverse_complement())
        config["rev_primer_rc"] = str(Seq(config["rev_primer"], IUPACAmbiguousDNA()).reverse_complement())
        
        shell("%(cutadapt)s -a forward=%(fwd_primer_rc)s -O {params.minovl} --trimmed-only -o primers/{wildcards.sample}_reverse_forward.fastq primers/{wildcards.sample}_reverse.fastq &>> {log}" % config)
        shell("%(cutadapt)s -a reverse=%(rev_primer_rc)s -O {params.minovl} --trimmed-only -o primers/{wildcards.sample}_forward_reverse.fastq primers/{wildcards.sample}_forward.fastq &>> {log}" % config)
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open("primers/%s_reverse_forward.fastq" % wildcards.sample), 
                                   "fastq"):
                newRec = rec.reverse_complement()
                newRec.id = "%s/rc" % rec.id
                newRec.description = rec.description.split(" ", 1)[1]
                out.write(newRec.format("fastq"))
        shell("cat primers/{wildcards.sample}_forward_reverse.fastq >> {output}")

rule fastq2fasta:
    input: "primers/{sample}_primer.fastq"
    output: "primers/{sample}_minLen.fasta"
    log: "logs/{sample}_minLen.log"
    params: minLen=2000
    run:
        with open(output[0], "w") as out, open(log[0], "w") as logFile:
            for read in SeqIO.parse(open(input[0]), "fastq"):
                if len(read) < params.minLen:
                    logFile.write("Removed %s, length %i < %i" % (read.id, len(read), params.minLen))
                else:
                    out.write(read.format("fasta"))

rule readNumbers_raw:
    input: raw=expand("raw/{sample}.fastq", sample=samples)
    output: "rawReadNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            #raw
            for rawFileName in input.raw:
                i=0
                with open(rawFileName) as rawFile:
                    iter = SeqIO.parse(rawFile, "fastq")
                    while True:
                        try:
                            next(iter)
                        except StopIteration:
                            break
                        i+=1
                sample = rawFileName.rsplit("/", 1)[-1].split(".")[0]
                out.write("raw\t%s\t%i\n" % (sample, i))
rule readNumbers_maxLen:
    input: length=expand("lenFilter/{sample}_rightLen.fastq", sample=samples)
    output: "maxLenReadNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            #max lenght filter
            for lenFileName in input.length:
                i=0
                with open(lenFileName) as lenFile:
                    iter = SeqIO.parse(lenFile, "fastq")
                    while True:
                        try:
                            next(iter)
                        except StopIteration:
                            break
                        i+=1
                sample = lenFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
                out.write("maxLenFilter\t%s\t%i\n" % (sample, i))

rule readNumbers_qual:
    input: qual=expand("qualFilter/{sample}_goodQual.fastq", sample=samples)
    output: "qualReadNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            #quality filter
            for qualFileName in input.qual:
                i=0
                with open(qualFileName) as qualFile:
                    iter = SeqIO.parse(qualFile, "fastq")
                    while True:
                        try:
                            next(iter)
                        except StopIteration:
                            break
                        i+=1
                sample = qualFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
                out.write("qualFilter\t%s\t%i\n" % (sample, i))

rule readNumbers_primer:
    input: primer=expand("primers/{sample}_primer.fastq", sample=samples)
    output: "primerReadNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            #primer filter
            for primerFileName in input.primer:
                i=0
                with open(primerFileName) as primerFile:
                    iter = SeqIO.parse(primerFile, "fastq")
                    while True:
                        try:
                            next(iter)
                        except StopIteration:
                            break
                        i+=1
                sample = primerFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
                out.write("primerFilter\t%s\t%i\n" % (sample, i))

rule readNumbers_minLen:
    input: minLen=expand("primers/{sample}_minLen.fasta", sample=samples)
    output: "minLenReadNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            #min length filter
            for minLenFileName in input.minLen:
                i=0
                with open(minLenFileName) as minLenFile:
                    iter = SeqIO.parse(minLenFile, "fasta")
                    while True:
                        try:
                            next(iter)
                        except StopIteration:
                            break
                        i+=1
                sample = minLenFileName.rsplit("/", 1)[-1].rsplit("_", 1)[0]
                out.write("minLenFilter\t%s\t%i\n" % (sample, i))

rule catReadNumber:
    input: "rawReadNumbers.tsv", "maxLenReadNumbers.tsv", "qualReadNumbers.tsv", "primerReadNumbers.tsv", "minLenReadNumbers.tsv"
    output: "readNumbers.tsv"
    shell: 
        "cat {input} > {output}"
        

rule plotReadNumber:
    input: "readNumbers.tsv"
    output: "readNumbers.pdf"
    run:
        R("""
        library(ggplot2)
        d=read.table("{input}")
        colnames(d) = c("stage", "sample", "number")
        d$stage=factor(d$stage, levels=c("raw", "maxLenFilter", "qualFilter", "primerFilter", "minLenFilter"))
        ggplot(d)+geom_bar(aes(x=sample, y=number, fill=stage), stat="identity", position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output}", width=16, height=10)
        """)

rule concat:
    input: expand("primers/{sample}_primer.fastq", sample=samples)
    output: fasta="primers/all_minLen.fasta", tsv="sampleInfo.tsv"
    log: "logs/all_minLen.log"
    params: minLen=2000
    run:
        with open(output.fasta, "w") as fasta, open(output.tsv, "w") as out, open(log[0], "w") as logFile:
            for inFile in input:
                for rec in SeqIO.parse(open(inFile), "fastq"):
                    if len(rec) < params.minLen:
                        logFile.write("Removed %s, length %i < %i" % (rec.id, len(rec), params.minLen))
                    else:
                        fasta.write(rec.format("fasta"))
                        sample = inFile.split("/")[-1].rsplit("_",1)[0]
                        out.write("%s\t%s\n" % (rec.id, sample))

rule cluster:
    input: "primers/{sample}_minLen.fasta"
    output: fasta="clusters/{sample}_cluster.fasta", clstr="clusters/{sample}_cluster.fasta.clstr"
    log: "logs/{sample}_clustering.log"
    threads: 3
    params: mem=16000
    shell:
        "%(cd-hit)s -i {input} -o {output.fasta} -g 1 -c 0.99 -r 0 -d 0 -M {params.mem} -T {threads} &> {log}" % config
        
rule clusterPerSample:
    input: cls="clusters/all_cluster.fasta.clstr", sample="sampleInfo.tsv"
    output: "clusters/all_cluster_persample.tsv"
    run:
        read2sample = {}
        for line in open(input.sample):
            read, sample = line.strip().split("\t")
            read2sample[read] = sample
        
        clusters = []
        for line in open(input.cls):
            if line[0] == ">":
                #new cluster
                clusters.append({})
            else:
                read = line.split("\t")[1].split(" ")[1][1:].strip(".")
                sample = read2sample[read]
                try:
                    clusters[-1][sample].append(read)
                except KeyError:
                    clusters[-1][sample] = [read]
        with open(output[0], "w") as out:
            for c, data in enumerate(clusters):
                for sample, seqs in data.items():
                    for read in seqs:
                        out.write("%s\t%s\tCluster %i\n" % (read, sample, c))

rule clusterSize:
    input: clsInfo="clusters/{sample}_cluster.fasta.clstr", clsFasta="clusters/{sample}_cluster.fasta"
    output: "clusters/{sample}_cluster.size.tsv"
    run:
        seqs = [r.id for r in SeqIO.parse(open(input.clsFasta), "fasta")]
        clusters = []
        for line in open(input.clsInfo):
            if line[0] == ">":
                #new cluster
                clusters.append(0)
            else:
                clusters[-1] += 1
        with open(output[0], "w") as out:
            for c, size in enumerate(clusters):
                out.write("%s\tCluster %i\t%s\t%i\n" % (wildcards.sample, c, seqs[c], size))

rule consensus:
    input: clusterTab="clusters/{sample}_cluster.fasta.clstr", clusterSeq="clusters/{sample}_cluster.fasta", reads="primers/{sample}_minLen.fasta"
    output: consensus="consensus/{sample}_consensus.fasta", reads="clusters/{sample}_cluster_reads.pic"
    log: aln="logs/{sample}_align.log", cons="logs/{sample}_consensus.log"
    params: minSize=3
    threads: 6
    run:
        minSize = params.minSize
        clusterReads = {}
        tReads = None
        repSeq = None
        for line in open(input.clusterTab):
            if line[0] == ">":
                #new cluster
                if not tReads is None:
                    assert not repSeq is None
                    clusterReads[repSeq] = tReads
                tReads = []
                repSeq = None
            else:
                #add to cluster
                redNum, seqInfo = line.strip().split("\t")
                lenStr, idStr, matchStr = seqInfo.split(" ", 2)
                readId = idStr[1:-3]
                tReads.append(readId)
                if matchStr == "*":
                    repSeq = readId
        #last cluster
        if not tReads is None:
            assert not repSeq is None
            clusterReads[repSeq] = tReads
        with open(output.reads, "wb") as pOut:
            pickle.dump(clusterReads, pOut)
        readRecs = {r.id: r for r in SeqIO.parse(open(input.reads), "fasta")}
        with open(output.consensus, "w") as out, open(log.cons, "w") as consLog:
            for clusterRec in SeqIO.parse(open(input.clusterSeq), "fasta"):
                size = len(clusterReads[clusterRec.id])
                if size < minSize:
                    continue
                baseFileName = "consensus/tmp_%s_%s" % (wildcards.sample, clusterRec.id.replace("/","_"))
                #write all reads in the cluster to temp file
                readFilePath = "%s_reads.fasta" % baseFileName
                with open(readFilePath, "w") as readFile:
                    for readId in clusterReads[clusterRec.id]:
                        readFile.write(readRecs[readId].format("fasta"))
                alignmentPath="%s_align.fasta" % baseFileName
                #run alignment
                shell("%(mafft)s" % config + " --ep 1 --thread {threads} %s > %s 2> {log.aln}" % (readFilePath, alignmentPath))
                #consensus
                consLog.write("======== %s\n" % clusterRec.id)
                cons = fastaConsensus(open(alignmentPath), log=consLog)
                out.write(">%s;size=%i;\n%s\n" % (clusterRec.id, size, cons))
#        shell("rm consensus/tmp_*")


rule prepChimeraRemoval:
    input: cls="clusters/{sample}_cluster.fasta", clsTab="clusters/{sample}_cluster.fasta.clstr"
    output: "clusters/{sample}_cluster_size.fasta"
    run:
        clusterSize = []
        rep=[]
        for line in open(input.clsTab):
            if line[0] == ">":
                #new cluster
                clusterSize.append(0)
            else:
                lenStr, nameStr, infoStr=line.strip().split("\t")[1].split(" ", 2)
                name = nameStr[1:-3]
                if infoStr == "*":
                    rep.append(name)
                clusterSize[-1] += 1
        sizeDict = dict(zip(rep, clusterSize))
        with open(output[0], "w") as out:
            for r, rec in enumerate(SeqIO.parse(open(input.cls), "fasta")):
                if sizeDict[rec.id] <2:
                    continue
                rec.id = "%s;size=%i;" % (rec.id, sizeDict[rec.id])
                rec.description = ""
                out.write(rec.format("fasta"))

rule removeChimera:
    input: seqs="consensus/{sample}_consensus.fasta"
    output: fasta="chimera/{sample}.nochimera.fasta", tsv="chimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_chimera.log"
    shell:
        "%(vsearch)s --uchime_denovo {input.seqs} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule removeChimeraNoCons:
    input: seqs="clusters/{sample}_cluster_size.fasta"
    output: fasta="chimeraNoCons/{sample}.nochimera.fasta", tsv="chimeraNoCons/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_chimeraNoCons.log"
    shell:
        "%(vsearch)s --uchime_denovo {input.seqs} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule nonChimeraReads:
    input: noChim="chimera/{sample}.nochimera.fasta", reads="primers/{sample}_minLen.fasta", cluster2read="clusters/{sample}_cluster_reads.pic"
    output: "chimera/{sample}.nochimeraReads.fasta"
    run:
        goodReads = {}
        with open(input.cluster2read, "rb") as pIn:
            clu2read = pickle.load(pIn)
        for rec in SeqIO.parse(open(input.noChim), "fasta"):
            try:
                readName = rec.id.split(";")[0]
                for readId in clu2read[readName]:
                    goodReads[readId] = None
            except KeyError:
                pass
        with open(output[0], "w") as out:
            for readRec in SeqIO.parse(input.reads, "fasta"):
                if readRec.id.split(";")[0] in goodReads:
                    out.write(readRec.format("fasta"))

rule itsx:
    input: "chimera/{sample}.nochimeraReads.fasta"
#    input: "chimera/{sample}.nochimera.fasta"
    output: "itsx/{sample}.SSU.fasta", "itsx/{sample}.ITS1.fasta", "itsx/{sample}.5_8S.fasta", "itsx/{sample}.ITS2.fasta", "itsx/{sample}.LSU.fasta", "itsx/{sample}.summary.txt", "itsx/{sample}.positions.txt", "itsx/{sample}.full.fasta"
    threads: 6
    log: "logs/{sample}_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 500 2> {log}" % config


rule otuCluster:
    input: "itsx/{sample}.full.fasta"
    output: fasta="otus/{sample}_otus.fasta", clstr="otus/{sample}_otus.fasta.clstr"
    log: "logs/{sample}_otuClustering.log"
    threads: 3
    params: mem=16000
    shell:
        "%(cd-hit)s -i {input} -o {output.fasta} -g 1 -c 0.97 -r 0 -d 0 -M {params.mem} -T {threads} &> {log}" % config

rule otuReads:
    input: clsInfo="otus/{sample}_otus.fasta.clstr"
    output: size="otus/{sample}_cluster.size.tsv", info="otus/{sample}_otuInfo.pic"
    run:
        clusterReads={}
        tReads = None
        repSeq = None
        for line in open(input.clsInfo):
            if line[0] == ">":
                #new cluster
                if not tReads is None:
                    assert not repSeq is None
                    clusterReads[repSeq] = tReads
                tReads = []
                repSeq = None
            else:
                #add to cluster
                redNum, seqInfo = line.strip().split("\t")
                lenStr, idStr, matchStr = seqInfo.split(" ", 2)
                readId = "%s/ITS" % idStr[1:-3].split("|")[0]
                tReads.append(readId)
                if matchStr == "*":
                    repSeq = readId
        #last cluster
        if not tReads is None:
            assert not repSeq is None
            clusterReads[repSeq] = tReads
        otuInf = {}
        with open(output.size, "w") as sOut:
            for otu, readList in clusterReads.items():
                for read in readList:
                    otuInf[read.rsplit("/", 1)[0]] = otu
                sOut.write("%s\t%i\n" % (otu, len(readList)))
        with open(output.info, "wb") as pOut:
            pickle.dump(otuInf, pOut)
        

rule transferOtus:
    input: pic="otus/{sample}_otuInfo.pic", fasta="itsx/{sample}.{marker}.fasta"
    output: "otus/{sample}_otus_{marker}.fasta"
    run:
        read2otu = pickle.load(open(input.pic, "rb"))
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input.fasta), "fasta"):
                readId = rec.id.split("|")[0]
                if readId in read2otu:
                    rec.id = "%s/%s" % (readId, wildcards.marker)
                    out.write(rec.format("fasta"))

####################################################################
# DBs

rule getUniteFile:
    output: "%(dbFolder)s/sh_general_release_dynamic_%(uniteVersion)s.fasta" % config
    shell:
        "cd %(dbFolder)s;" \
        "wget https://unite.ut.ee/sh_files/sh_general_release_%(uniteVersion)s.zip;" \
        "unzip sh_general_release_%(uniteVersion)s.zip;" \
        "rm sh_general_release_%(uniteVersion)s.zip" % config

rule createUniteTax:
    input: "%(dbFolder)s/sh_general_release_dynamic_%(uniteVersion)s.fasta" % config
    output: fasta="%(dbFolder)s/UNITE_%(uniteVersion)s.fasta" % config, tax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
    run:
        with open(output.tax, "w") as tOut, open(output.fasta, "w") as fOut:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                arr = rec.id.split("|")
                sh = arr[2]
                tax = arr[-1]
                tOut.write("%s\t%s\n" % (sh, tax))
                rec.id = sh
                fOut.write(rec.format("fasta"))

rule creatUniteIndex:
    input: "%(dbFolder)s/UNITE_%(uniteVersion)s.fasta" % config
    output: touch("%(dbFolder)s/UNITE_%(uniteVersion)s.fasta.lambdaIndexCreated" % config)
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config
        
rule getSilva_main:
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz" % config
    shell:
        "cd %(dbFolder)s;" \
        "wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz" % config

rule getSilva_md5:
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz.md5" % config
    shell: 
        "cd %(dbFolder)s;" \
        "wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz.md5" % config
    
rule getSilva_test:
    input: gz="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz" % config, md5="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz.md5" % config
    output: touch("%(dbFolder)s/silva_{marker}dl_good" % config)
    shell: 
        "cd %(dbFolder)s;" \
        "md5sum -c SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz.md5" % config

rule unpackSilva:
    input: gz="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz" % config, good="%(dbFolder)s/silva_{marker}dl_good" % config
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta" % config
    shell:
        "cd %(dbFolder)s;" \
        "gunzip SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz; " \
        "touch SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta" % config

rule createSlivaTax:
    input: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta" % config
    output: tax="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}_tax.tsv" % config
    run:
        with open(output.tax, "w") as tOut:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                tax = rec.description.split(" ", 1)[1].replace(" ", "_")
                tOut.write("%s\t%s\n" % (rec.id, tax))

rule creatSilvaIndex:
    input: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta" % config
    output: touch("%(dbFolder)s/silva_{marker}_lambdaIndexCreated" % config)
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config

####################################################################


rule alignToUnite:
    input: clu="otus/{sample}_otus.fasta", db="%(dbFolder)s/UNITE_%(uniteVersion)s.fasta" % config, dbFlag="%(dbFolder)s/UNITE_%(uniteVersion)s.fasta.lambdaIndexCreated" % config
    output: "lambda/{sample}.otu_vs_UNITE.m8"
    log: "logs/{sample}_otu_lambda.log"
    threads: 3
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -d {input.db} -o {output} -p blastn -t {threads} &> {log}" % config

rule classifyITS:
    input: lam="lambda/{sample}.otu_vs_UNITE.m8", clu="otus/{sample}_otus.fasta", tax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
    output: "taxonomy/{sample}_otu_ITS.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/{sample}_otuClass.log"
    run:
        taxDict = {}
        for line in open(input.tax):
            sh, tax = line.strip().split("\t")
            taxDict[sh] = tax
        logOut = open(log[0], "w")
        classifi = {}
        seqLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
#        for rec in SeqIO.parse(open(input.clu), "fasta"):
#            seqNr += 1
#            classifi[rec.id] = []
#            seqLength[rec.id] = len(rec)
        for line in open(input.lam, encoding="latin-1"):
            total +=1
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
            readId = qseqid.split("|")[0]
            if float(evalue) > params.maxE:
                evalueFilter += 1
                continue
            if float(pident) < params.minIdent:
                identFilter +=1
                continue
#            if float(length)/seqLength[readId]*100 < params.minCov:
#                covFilter += 1
#                continue
            linStr = taxDict[sseqid]
            try:
                classifi[readId].append((linStr, float(bitscore)))
            except KeyError:
                classifi[readId] = [(linStr, float(bitscore))]
        logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
        logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, params.maxE))
        logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, params.minIdent))
        logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, params.minCov))
        topPerc = params.topPerc/100.0
        with open(output[0], "w") as out:
            for key, hits in classifi.items():
                if not hits:
                    out.write("%s\tunknown\n" % (key))
                else:
                    sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                    cutoff = 0
                    while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                        cutoff += 1
                    lineage = lca([hit[0] for hit in sortedHits[:cutoff]], params.stringency)
                    out.write("%s\t%s\n" % (key, lineage))
        try:
            logOut.close()
        except:
            pass

rule alignToSliva:
    input: clu="otus/{sample}_otus_{marker}.fasta", db="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta" % config, dbFlag="%(dbFolder)s/silva_{marker}_lambdaIndexCreated" % config
    output: "lambda/{sample}.{marker}_vs_SILVA.m8"
    log: "logs/{sample}_otu{marker}_lambda.log"
    threads: 3
    shell:
        "%(lambdaFolder)s/lambda -q {input.clu} -d {input.db} -o {output} -p blastn -t {threads} &> {log}" % config

rule classifySilva:
    input: lam="lambda/{sample}.{marker}_vs_SILVA.m8", clu="otus/{sample}_otus.fasta", tax="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}_tax.tsv" % config
    output: "taxonomy/{sample}_otu_{marker}.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/{sample}_{marker}_otuClass.log"
    run:
        taxDict={}
        for line in open(input.tax):
            rId, tax = line.strip().split("\t")
            taxDict[rId] = tax
        logOut = open(log[0], "w")
        classifi = {}
        seqLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
#        for rec in SeqIO.parse(open(input.clu), "fasta"):
#            seqNr += 1
#            classifi[rec.id] = []
#            seqLength[rec.id] = len(rec)
        for line in open(input.lam, encoding="latin-1"):
            total +=1
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
            readId = qseqid
            if float(evalue) > params.maxE:
                evalueFilter += 1
                continue
            if float(pident) < params.minIdent:
                identFilter +=1
                continue
#            if float(length)/seqLength[readId]*100 < params.minCov:
#                covFilter += 1
#                continue
            linStr = taxDict[sseqid]
            try:
                classifi[qseqid].append((linStr, float(bitscore)))
            except KeyError:
                classifi[qseqid]= [(linStr, float(bitscore))]
        logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
        logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, params.maxE))
        logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, params.minIdent))
        logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, params.minCov))
        topPerc = params.topPerc/100.0
        with open(output[0], "w") as out:
            for key, hits in classifi.items():
                if not hits:
                    out.write("%s\tunknown\n" % (key))
                else:
                    sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                    cutoff = 0
                    while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                        cutoff += 1
                    lineage = lca([hit[0] for hit in sortedHits[:cutoff]], params.stringency)
                    out.write("%s\t%s\n" % (key, lineage))
        try:
            logOut.close()
        except:
            pass


rule compareCls:
    input: ssu="taxonomy/{sample}_otu_SSU.class.tsv", its="taxonomy/{sample}_otu_ITS.class.tsv", lsu="taxonomy/{sample}_otu_LSU.class.tsv", read2otu="otus/{sample}_otuInfo.pic"
    output: "taxonomy/{sample}_comb.class.tsv"
    params: stringency=.90
    run:
        read2otu = pickle.load(open(input.read2otu, "rb"))
        otu2read = {}
        for read, otu in read2otu.items():
            try:
                otu2read[otu].append(read)
            except KeyError:
                otu2read[otu] = [read]
        ssu = {}
        for line in open(input.ssu):
            ssuId, ssuTax = line.strip().split("\t")
            ssu[ssuId] = ssuTax
        lsu = {}
        for line in open(input.lsu):
            lsuId, lsuTax = line.strip().split("\t")
            lsu[lsuId] = lsuTax
        with open(output[0], "w") as out:
            for line in open(input.its):
                itsId, itsTax = line.strip().split("\t")
                otuReads = otu2read["%s/ITS" % itsId]
                lsuTax = lca([lsu["%s/LSU" %r] for r in otuReads], params.stringency)
                ssuTax = lca([ssu["%s/SSU" %r] for r in otuReads], params.stringency)
                out.write("%s\t%s\t%s\n" % (ssuTax, itsTax, lsuTax))

def lca(lineageStrings, stringency=1.0, 
        unidentified=["unidentified", "unclassified", "unknown"],
        ignoreIncertaeSedis=True):
    lineage = []
    mLineages = []
    #remove bootstrap values ("(100)", "(75)", etc.) if any
    for mLin in [l.strip(";").split(";") for l in lineageStrings]:
        mLineages.append([])
        for entry in mLin:
             mLineages[-1].append(entry.split("(")[0])
    i=0
    maxLinLen = max([len(m) for m in mLineages])
    active = [True]*len(mLineages)
    for i in range(maxLinLen):
        total = 0.0
        counts = {}
        for m, memberLin in enumerate(mLineages):
            if not active[m]:
                continue #ignore lineages that were deactivated further up in the tree
            if len(memberLin) <= i:
                active[m] = False
                continue #ignore lineages that are not this long
            name = memberLin[i].split("__")[-1]
            if name in unidentified:
                continue # ignoring unidentified entrys
            if ignoreIncertaeSedis and name.startswith("Incertae"):
                continue # ignoring Incertae sedis entries.
                         # NOTE: this will mean lineages end at the first Incerta sedis
            total += 1
            try:
                counts[memberLin[i]] += 1
            except KeyError:
                counts[memberLin[i]] = 1
        if not counts:
            #no valid lineage entrys found in this level
            break
        most=sorted(counts.items(), key=lambda x: x[1], reverse=True)[0]
        different = total - most[1]
        #accept the lineage entry if its proportion of all (valid) classifications higher than stringency setting
        if different/total <= (1.0-stringency):
            lineage.append(most[0]) #add the most apearing entry to the new lineage
            #deactivate all lineages that were different at this level
            for m, memberLin in enumerate(mLineages):
                if active[m] and memberLin[i] != most[0]:
                    active[m] = False
        else:
            break
    if len(lineage) == 0:
        lineage = ["unknown"]
    return ";".join(lineage)


#####################################################################

#rule cluster2:
#    input: "consensus/{sample}_consensus.fasta"
#    output: fasta="clusters2/{sample}_cluster2.fasta", clstr="clusters2/{sample}_cluster2.fasta.clstr"
#    log: "logs/{sample}_clustering2.log"
#    threads: 3
#    shell:
#        "%(cd-hit)s -i {input} -o {output.fasta} -g 1 -c 0.99 -r 0 -d 0 -T {threads} &> {log}" % config
#        
#rule cluster2Sizes:
#    input: clsInfo="clusters2/{sample}_cluster2.fasta.clstr"
#    output: "clusters2/{sample}_cluster2.size.tsv"
#    run:
#        clusters = []
#        rep=[]
#        for line in open(input.clsInfo):
#            if line[0] == ">":
#                #new cluster
#                clusters.append(0)
#            else:
#                lenStr, nameStr, infoStr=line.strip().split("\t")[1].split(" ", 2)
#                name, sizeStr, _ = nameStr[1:].split(";")
#                size = int(sizeStr.split("=")[1])
#                if infoStr == "*":
#                    rep.append(name)
#                clusters[-1] += size
#        with open(output[0], "w") as out:
#            for c, size in enumerate(clusters):
#                out.write("%s\tCluster %i\t%s\t%i\n" % (wildcards.sample, c, rep[c], size))

#rule consensus2:
#    input: clusterTab="clusters2/{sample}_cluster2.fasta.clstr", clusterSeq="clusters2/{sample}_cluster2.fasta", reads="consensus/{sample}_consensus.fasta"
#    output: consensus="consensus2/{sample}_consensus2.fasta"
#    log: cons="logs/{sample}_consensus2.log", aln="logs/{sample}_align2.log"
#    params: minSize=2
#    threads: 6
#    run:



def fastaConsensus(inFasta, warn=0.75, log=sys.stderr):
    prof = []
    for rec in SeqIO.parse(inFasta, "fasta"):
        for c, char in enumerate(rec.seq):
            try:
                prof[c][char.upper()] += 1
            except IndexError:
                prof.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
                prof[c][char.upper()] += 1
    cons = []
    for i, col in enumerate(prof):
        cov = float(sum(col.values()))
        most = sorted(list(col.items()), key=lambda x: x[1], reverse=True)[0]
        if most[1]/cov < warn:
            log.write("At position %i the most common base %s only has a frequency of %f\n" % (i, most[0], most[1]/cov))
        if most[0] != "-":
            cons.append(most[0])
    return "".join(cons)

#def samConsensus(inSam, warn=0.75, log=sys.stderr):
#    prof = []
#    for line in inSam:
#        if line[0] == "@":
#            continue
#        aSeq = []
#        readId, flag, refId, pos, mapQ, cigar, pnext, rnext, alen, seq, seqQ, tags = line.strip().split("\t", 11)
#        aSeq.append("-"*(int(pos)-1))
#        cur=0
#        for entry in re.findall("[0-9]+[MIDNSHP=X]", cigar):
#            l = int(entry[:-1])
#            if entry[-1] == "M":
#                aSeq.append(seq[cur:cur+l])
#                cur += l
#            elif entry[-1] == "I":
#                pass
#                cur += l
#            elif entry[-1] == "D":
#                aSeq.append("-"*l)
#            else:
#                raise NotImplementedError("Cigar entry %s not implemented yet" % entry[-1])
##        log.write(">%s\n%s\n" % (readId, "".join(aSeq)))
#        for c,char in enumerate("".join(aSeq)):
#            try:
#                prof[c][char] += 1
#            except IndexError:
#                prof.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
#                prof[c][char] += 1
#    cons = []
#    for i, col in enumerate(prof):
#        cov = float(sum(col.values()))
#        most = sorted(list(col.items()), key=lambda x: x[1], reverse=True)[0]
#        if most[1]/cov < warn:
#            log.write("At position %i the most common base %s only has a frequency of %f\n" % (i, most[0], most[1]/cov))
#        if most[0] != "-":
#            cons.append(most[0])
#    return "".join(cons)


