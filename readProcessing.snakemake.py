rule unpack:
    input: "%(inFolder)s/8_libs_Mar17/Ampl.Lib{cellNr}.SC1+2_barcoded-fastqs.tgz" % config
    output: dynamic("%(inFolder)s/Lib{cellNr}_1+2/{barcode}_Forward--{barcode}_Forward.fastq" % config)
    shell:
        "mkdir -p %(inFolder)s/Lib{wildcards.cellNr}; tar -xzf {input} -C %(inFolder)s/Lib{wildcards.cellNr}_1+2; touch %(inFolder)s/Lib{wildcards.cellNr}_1+2/*.fastq" % config

rule unpackThird:
    input: "%(inFolder)s/8_libs_Mar17/Ampl.Lib{cellNr}.SC3_barcoded-fastqs.tgz" % config
    output: dynamic("%(inFolder)s/Lib{cellNr}_3/{barcode}_Forward--{barcode}_Forward.fastq" % config)
    shell:
        "mkdir -p %(inFolder)s/Lib{wildcards.cellNr}; tar -xzf {input} -C %(inFolder)s/Lib{wildcards.cellNr}_3; touch %(inFolder)s/Lib{wildcards.cellNr}_3/*.fastq" % config

rule concatRawfiles:
    input: "%(inFolder)s/Lib{cellNr}_1+2/{barcode,\d+}_Forward--{barcode,\d+}_Forward.fastq" % config, "%(inFolder)s/Lib{cellNr}_3/{barcode,\d+}_Forward--{barcode,\d+}_Forward.fastq" % config
    output: "%(inFolder)s/Lib{cellNr}-{barcode,\d+}.fastq" % config
    shell:
        "cat {input} > {output}"

rule fastqc:
    input: "%(inFolder)s/{sample}.fastq" % config
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
    output: right="lenFilter/{sample}_rightLen.fastq", long="lenFilter/{sample}_tooLong.fastq", short="lenFilter/{sample}_tooShort.fastq"
    log: "logs/{sample}_lenFilter.log"
    params: maxLen = 6500, minLen = 3000
    run:
        nrLong = 0
        nrShort = 0
        with open(output.right, "w") as out, open(output.long, "w") as tooLong, open(output.short, "w") as tooShort:
            for read in SeqIO.parse(open(input[0]), "fastq"):
                if len(read) <= params.minLen:
                    tooShort.write(read.format("fastq"))
                    nrShort += 1
                elif len(read) <= params.maxLen:
                    out.write(read.format("fastq"))
                else:
                    tooLong.write(read.format("fastq"))
                    nrLong += 1
        with open(log[0], "w") as logFile:
            logFile.write("%i reads were removed because they were longer than %i\n" % (nrLong, params.maxLen))
            logFile.write("%i reads were removed because they were shorter than %i\n" % (nrLong, params.maxLen))

rule qualityFilter:
    input: fastq="lenFilter/{sample}_rightLen.fastq"
    output: good="qualFilter/{sample}_goodQual.fastq", bad="qualFilter/{sample}_badQual.fastq", info="qualFilter/{sample}_qualInfo.tsv"
    params: minQual=0.996
    log: "logs/{sample}_qualityFilter.log"
    run:
        removed = 0
        
        with open(output.good, "w") as out, open(output.bad, "w") as trash, open(output.info, "w") as info:
            for read in SeqIO.parse(open(input.fastq), "fastq"):
                try:
                    tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"]]) / len(read)
                except Exception:
                    print(read.id)
                    raise
                info.write("%s\t%f\n" % (read.id, 1-tError))
                if (1-tError) < params.minQual:
                    removed += 1
                    trash.write(read.format("fastq"))
                else:
                    out.write(read.format("fastq"))
        open(log[0], "w").write("%s: %i reads removed because quality < %f\n" % (wildcards.sample, removed, params.minQual))

rule qualityVsLength:
    input: fasta="lenFilter/{sample}_rightLen.fastq", qual="qualFilter/{sample}_qualInfo.tsv"
    output: tsv="qualFilter/{sample}_LenVsQual.tsv"
    run:
        qual = {}
        for line in open(input.qual):
            rId, q = line.strip().split("\t")
            qual[rId] = float(q)
        with open(output.tsv, "w") as out:
            for rec in SeqIO.parse(open(input.fasta), "fastq"):
                out.write("%s\t%s\t%i\t%f\n" % (wildcards.sample, rec.id, len(rec), qual[rec.id]))

rule concatQualVsLength:
    input: expand("qualFilter/{sample}_LenVsQual.tsv", sample=samples)
    output: "qualFilter/lenVsQual.tsv"
    shell:
        "cat {input} > {output}"

rule plotQualVsLength:
    input: "qualFilter/lenVsQual.tsv"
    output: "lenVsQual.png"
    run:
        R("""
        library(ggplot2)
        d=read.table("{input}", header=F)
        colnames(d) = c("sample", "read", "length", "qual")
        m=lm(qual~length+sample, d)
        
        ggplot(d) + geom_point(aes(length, qual, color=sample), alpha=0.5) + geom_smooth(aes(length, qual), method = "lm", linetype = 2)
        ggsave("{output}", width=16, height=10)
        """)
        

rule windowQualFilter:
    input: fastq="qualFilter/{sample}_goodQual.fastq"
    output: good="windowQualFilter/{sample}_goodQual.fastq", stat="windowQualFilter/{sample}_stat.tsv"
    params: winSize=8, minQual=0.9
    log: "logs/{sample}_winQualityFilter.log"
    run:
        removed = 0
        with open(output.good, "w") as out, open(output.stat, "w") as statOut:
            for read in SeqIO.parse(open(input.fastq), "fastq"):
                anyRemoved = False
                for i in range(len(read)-params.winSize):
                    tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"][i:i+params.winSize] ]) / params.winSize
                    tQual = 1 - tError
                    tRemoved=False
                    if (tQual) < params.minQual:
                        removed += 1
                        tRemoved = True
                    kmer=str(read.seq[i:i+params.winSize])
                    hp = homopoly(kmer)
                    statOut.write("%s\t%s\t%i\t%i\t%f\t%s\t%s\t%i\t%i\t%i\t%i\t%i\n" % (wildcards.sample, read.id, i, len(read), tQual, tRemoved, kmer, kmer.count("A"), kmer.count("C"), kmer.count("G"), kmer.count("T"), hp))
                if not anyRemoved:
                    out.write(read.format("fastq"))
            with open(log[0], "w") as logFile:
                logFile.write("%s: %i reads removed because in a window of size %i quality droped below %f\n" % (wildcards.sample, removed, params.winSize, params.minQual))

rule catWindowQual:
    input: expand("windowQualFilter/Lib{libNr}-0075_stat.tsv", libNr=range(1,9)), expand("windowQualFilter/Lib{libNr}-0034_stat.tsv", libNr=[1,2,3,5,6,7,8]), "windowQualFilter/Lib4-0018_stat.tsv"
    output: "windowQualFilter/envStats.tsv"
    shell:
        "cat {input} > {output}"

rule plotWindowQual:
    input: "windowQualFilter/envStats.tsv"
    

rule filterPrimer:
    input: "windowQualFilter/{sample}_goodQual.fastq"
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
    output: "primers/{sample}_primer.fasta"
    run:
        with open(output[0], "w") as out:
            for read in SeqIO.parse(open(input[0]), "fastq"):
                out.write(read.format("fasta"))

rule readNumbers_raw:
    input: raw=expand("raw/{sample}.fastq", sample=samples)
    output: "readNumbers/rawReadNumbers.tsv"
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
    output: "readNumbers/lenReadNumbers.tsv"
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
                out.write("lenFilter\t%s\t%i\n" % (sample, i))

rule readNumbers_qual:
    input: qual=expand("qualFilter/{sample}_goodQual.fastq", sample=samples)
    output: "readNumbers/qualReadNumbers.tsv"
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

rule readNumbers_winQual:
    input: qual=expand("windowQualFilter/{sample}_goodQual.fastq", sample=samples)
    output: "readNumbers/winQualReadNumbers.tsv"
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
                out.write("winQualFilter\t%s\t%i\n" % (sample, i))

rule readNumbers_primer:
    input: primer=expand("primers/{sample}_primer.fastq", sample=samples)
    output: "readNumbers/primerReadNumbers.tsv"
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

rule catReadNumber:
    input: "readNumbers/rawReadNumbers.tsv", "readNumbers/lenReadNumbers.tsv", "readNumbers/qualReadNumbers.tsv", "readNumbers/winQualReadNumbers.tsv", "readNumbers/primerReadNumbers.tsv"
    output: "readNumbers/readNumbers.tsv"
    shell: 
        "cat {input} > {output}"
        

rule plotReadNumber:
    input: "readNumbers/readNumbers.tsv"
    output: "readNumbers.pdf"
    run:
        R("""
        library(ggplot2)
        d=read.table("{input}")
        colnames(d) = c("stage", "sample", "number")
        d$stage=factor(d$stage, levels=c("raw", "lenFilter", "qualFilter", "winQualFilter", "primerFilter"))
        ggplot(d)+geom_bar(aes(x=sample, y=number, fill=stage), stat="identity", position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave("{output}", width=16, height=10)
        """)
        
rule prepPrecluster:
    input: fasta="primers/{sample}_primer.fasta", qual="qualFilter/{sample}_qualInfo.tsv"
    output: "preclusters/{sample}_cluInput.fasta"
    run:
        qual = {}
        for line in open(input.qual):
            rId, tQual = line.strip().split("\t")
            qual[rId] = tQual
        readList = []
        for read in SeqIO.parse(open(input.fasta), "fasta"):
            if read.id.endswith("rc"):
                readId = read.id.rsplit("/", 1)[0]
            else:
                readId = read.id
            readList.append((read, qual[readId]))
        with open(output[0], "w") as out:
            for read, qual in sorted(readList, key=lambda x: x[1], reverse=True):
                read.description=str(qual)
                out.write(read.format("fasta"))

rule preCluster:
    input: "preclusters/{sample}_cluInput.fasta"
    output: uc="preclusters/{sample}.uc.txt", cons="consensus/{sample}_consensus.fasta"
    log: "logs/{sample}_pre-cluster.log"
    threads: 6
    shell:
        "%(vsearch)s --usersort --cluster_smallmem {input} --relabel {wildcards.sample}_precluster --sizeout --iddef 0 --id 0.99 --minsl 0.9 --consout {output.cons} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

rule preClusterInfo:
    input: clsInfo="preclusters/{sample}.uc.txt"
    output: info="preclusters/{sample}_cluInfo.tsv", size="preclusters/{sample}_cluster.size.tsv"
    run:
        with open(output.info, "w") as out, open(output.size, "w") as sOut:
            for line in open(input.clsInfo):
                if line[0] == "C":
                    rType, cNr, size, _ = line.split("\t", 3)
                    sOut.write("precluster%i\t%s\n" % (int(cNr)+1, size))
                elif line[0] in "SH":
                    arr = line.strip().split("\t")
                    seq = arr[-2]
                    cluNr = int(arr[1])
                    out.write("%s_precluster%i\t%s\n" % (wildcards.sample, cluNr+1, seq))
                else:
                    raise ValueError("Unknown record type: %s" % line[0])

rule removeChimera:
    input: seqs="consensus/{sample}_consensus.fasta"
    output: fasta="chimera/{sample}.nochimera.fasta", tsv="chimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_chimera.log"
    shell:
        "%(vsearch)s --uchime_denovo {input.seqs} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule itsx:
    input: "chimera/{sample}.nochimera.fasta"
    output: "itsx/{sample}.SSU.fasta", "itsx/{sample}.ITS1.fasta", "itsx/{sample}.5_8S.fasta", "itsx/{sample}.ITS2.fasta", "itsx/{sample}.LSU.fasta", "itsx/{sample}.summary.txt", "itsx/{sample}.positions.txt", "itsx/{sample}.full.fasta"
    threads: 6
    log: "logs/{sample}_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 500 -E 1e-4 2> {log}" % config

def sampleMappingInput(wildcards):
    if wildcards.sampleSet == "stechlin":
        return expand("itsx/{sampleSet}.full.fasta", sampleSet=stechlin)
    else:
        return expand("itsx/{sampleSet}.full.fasta", sampleSet=samples)

rule getSampleMapping:
    input: sampleMappingInput
    output: sample="{sampleSet}_preClu2sample.pic"
    run:
        preClu2sample = {}
        for inputFile in input:
            sample = inputFile.rsplit("/", 1)[-1].split(".", 1)[0]
            for rec in SeqIO.parse(open(inputFile), "fasta"):
                preClu2sample[rec.id] = sample
        pickle.dump(preClu2sample, open(output.sample, "wb"))

def homopoly(kmer):
    last=None
    rv=1
    h=1
    for j in range(len(kmer)):
        if kmer[j] == last:
            h +=1
        else:
            rv=max(rv, h)
            h=1
            last = kmer[j]
    return(max(rv, h))
