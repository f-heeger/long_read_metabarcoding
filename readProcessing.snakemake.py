
rule getData:
    """Download read data from the SRA"""
    output: "%(inFolder)s/{sample}.fastq" % config
    run:
        sId = sampleInfo["sraID"][wildcards.sample]
        shell("%(fastq-dump)s" % config + " %s -Z > {output}" % sId)

rule fastqc:
    """Run fastqc an all samples"""
    input: "%(inFolder)s/{sample}.fastq" % config
    output: "QC/{sample}_fastqc.html"
    threads: 6
    shell: "%(fastqc)s -o QC -t {threads} {input}" % config

rule multiqc:
    """Combine QC reports into one with multiqc"""
    input: expand("QC/{sample}_fastqc.html", sample=allSamples)
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    shell:
        "%(multiqc)s -f --interactive -o QC QC/*_fastqc.zip" % config

rule lengthFilter:
    """Filter reads by minimum and maximum length"""
    input: fastq="%(inFolder)s/{sample}.fastq" % config
    output: right="lenFilter/{sample}_rightLen.fastq", long="lenFilter/{sample}_tooLong.fastq", short="lenFilter/{sample}_tooShort.fastq"
    log: "logs/{sample}_lenFilter.log"
    run:
        maxLen = config["maxReadLen"]
        minLen = config["minReadLen"]
        nrLong = 0
        nrShort = 0
        with open(output.right, "w") as out, open(output.long, "w") as tooLong, open(output.short, "w") as tooShort:
            for read in SeqIO.parse(open(input[0]), "fastq"):
                if len(read) <= minLen:
                    tooShort.write(read.format("fastq"))
                    nrShort += 1
                elif len(read) <= maxLen:
                    out.write(read.format("fastq"))
                else:
                    tooLong.write(read.format("fastq"))
                    nrLong += 1
        with open(log[0], "w") as logFile:
            logFile.write("%i reads were removed because they were longer than %i\n" % (nrLong, maxLen))
            logFile.write("%i reads were removed because they were shorter than %i\n" % (nrLong, maxLen))

rule qualityFilter:
    """Filter reads by minimum mean quality"""
    input: fastq="lenFilter/{sample}_rightLen.fastq"
    output: good="qualFilter/{sample}_goodQual.fastq", bad="qualFilter/{sample}_badQual.fastq", info="qualFilter/{sample}_qualInfo.tsv"
    log: "logs/{sample}_qualityFilter.log"
    run:
        minQual = config["minReadQual"]
        removed = 0
        
        with open(output.good, "w") as out, open(output.bad, "w") as trash, open(output.info, "w") as info:
            for read in SeqIO.parse(open(input.fastq), "fastq"):
                try:
                    tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"]]) / len(read)
                except Exception:
                    print(read.id)
                    raise
                info.write("%s\t%f\n" % (read.id, 1-tError))
                if (1-tError) < minQual:
                    removed += 1
                    trash.write(read.format("fastq"))
                else:
                    out.write(read.format("fastq"))
        open(log[0], "w").write("%s: %i reads removed because quality < %f\n" % (wildcards.sample, removed, minQual))

rule qualityVsLength:
    """Generate data for plotting length vs quality"""
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
    """Concatenate quality vs length data for all samples"""
    input: expand("qualFilter/{sample}_LenVsQual.tsv", sample=samples)
    output: "qualFilter/lenVsQual.tsv"
    shell:
        "cat {input} > {output}"

rule plotQualVsLength:
    """Plot quality vs length"""
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
    """Filter by sliding window mean quality"""
    input: fastq="qualFilter/{sample}_goodQual.fastq"
    output: good="windowQualFilter/{sample}_goodQual.fastq", stat="windowQualFilter/{sample}_stat.tsv"
    log: "logs/{sample}_winQualityFilter.log"
    run:
        winSize = config["qualityWindowSize"]
        minQual = config["windowMinQuality"]
        removed = 0
        with open(output.good, "w") as out, open(output.stat, "w") as statOut:
            for read in SeqIO.parse(open(input.fastq), "fastq"):
                anyRemoved = False
                for i in range(len(read)-winSize):
                    tError = sum([10.0**(float(-q)/10.0) for q in read.letter_annotations["phred_quality"][i:i+winSize] ]) / winSize
                    tQual = 1 - tError
                    tRemoved=False
                    if (tQual) < minQual:
                        tRemoved = True
                        anyRemoved = True
                    kmer=str(read.seq[i:i+winSize])
                    hp = homopoly(kmer)
                    statOut.write("%s\t%s\t%i\t%i\t%f\t%s\t%s\t%i\t%i\t%i\t%i\t%i\n" % (wildcards.sample, read.id, i, len(read), tQual, tRemoved, kmer, kmer.count("A"), kmer.count("C"), kmer.count("G"), kmer.count("T"), hp))
                if not anyRemoved:
                    out.write(read.format("fastq"))
                else:
                    removed += 1
            with open(log[0], "w") as logFile:
                logFile.write("%s: %i reads removed because in a window of size %i quality droped below %f\n" % (wildcards.sample, removed, winSize, minQual))

rule catWindowQual:
    """concatenate window quality data for all environmental samples"""
    input: expand("windowQualFilter/{sample}_stat.tsv", sample=samples)
    output: "windowQualFilter/envStats.tsv"
    shell:
        "cat {input} > {output}"

rule filterPrimer:
    """Filter sequences by occurence of primer sequences and cut primer sequences"""
    input: "windowQualFilter/{sample}_goodQual.fastq"
    output: fastq="primers/{sample}_primer.fastq"
    log: "logs/{sample}_primer.log"
    threads: 3
    run:
        shell("%(cutadapt)s -g forward=%(fwd_primer)s -g reverse=%(rev_primer)s -O %(primerMinOverlap)i --trimmed-only -o primers/{wildcards.sample}_{{name}}.fastq {input} &> {log}" % config)
        config["fwd_primer_rc"] = str(Seq(config["fwd_primer"], IUPACAmbiguousDNA()).reverse_complement())
        config["rev_primer_rc"] = str(Seq(config["rev_primer"], IUPACAmbiguousDNA()).reverse_complement())
        
        shell("%(cutadapt)s -a forward=%(fwd_primer_rc)s -O %(primerMinOverlap)i --trimmed-only -o primers/{wildcards.sample}_reverse_forward.fastq primers/{wildcards.sample}_reverse.fastq &>> {log}" % config)
        shell("%(cutadapt)s -a reverse=%(rev_primer_rc)s -O %(primerMinOverlap)i --trimmed-only -o primers/{wildcards.sample}_forward_reverse.fastq primers/{wildcards.sample}_forward.fastq &>> {log}" % config)
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open("primers/%s_reverse_forward.fastq" % wildcards.sample), 
                                   "fastq"):
                newRec = rec.reverse_complement()
                newRec.id = "%s/rc" % rec.id
                newRec.description = rec.description.split(" ", 1)[1]
                out.write(newRec.format("fastq"))
        shell("cat primers/{wildcards.sample}_forward_reverse.fastq >> {output}")

rule fastq2fasta:
    """Convert reads from fastq to fasta"""
    input: "primers/{sample}_primer.fastq"
    output: "primers/{sample}_primer.fasta"
    run:
        with open(output[0], "w") as out:
            for read in SeqIO.parse(open(input[0]), "fastq"):
                out.write(read.format("fasta"))
        shell("sleep 1m; touch {output}")

rule readNumbers_raw:
    """Count raw reads"""
    input: raw=expand("%(inFolder)s/{sample}.fastq" % config, sample=allSamples)
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
    """Count reads after length filter"""
    input: length=expand("lenFilter/{sample}_rightLen.fastq", sample=allSamples)
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
    """Count reads after average quality filter"""
    input: qual=expand("qualFilter/{sample}_goodQual.fastq", sample=allSamples)
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
    """Countr reads after sliding window quality filter"""
    input: qual=expand("windowQualFilter/{sample}_goodQual.fastq", sample=allSamples)
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
    """Count reads after primer filtering"""
    input: primer=expand("primers/{sample}_primer.fastq", sample=allSamples)
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
    """Concatenate all read numbers for plotting"""
    input: "readNumbers/rawReadNumbers.tsv", "readNumbers/lenReadNumbers.tsv", "readNumbers/qualReadNumbers.tsv", "readNumbers/winQualReadNumbers.tsv", "readNumbers/primerReadNumbers.tsv"
    output: "readNumbers/readNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            for inFile in input:
                for line in open(inFile):
                    stage, sample, number = line.strip().split("\t")
                    lib, bc = sample.split("-")
                    out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (stage, sample, sampleName[sample], number, lib, bc))

rule plotReadNumber:
    """Plot read number after each filtering step"""
    input: "readNumbers/readNumbers.tsv"
    output: graph="readNumbers.svg", tab="readNumbersTab.tsv"
    run:
        R("""
        library(ggplot2)
        library(reshape2)

        collectData <- function (set, setname, a_total, a_rel) {{
            set$prop=NA

            for (tSample in unique(set$sample)) {{
                for (stage in c("lenFilter", "qualFilter", "winQualFilter", "primerFilter")) {{
                    set[set$sample==tSample & set$stage==stage,]$prop = set[set$sample==tSample & set$stage==stage,]$number/set[set$sample==tSample & set$stage=="raw",]$number
                }}
                set[set$sample==tSample & set$stage=="raw",]$prop = 1.0
            }}

            new = aggregate(prop~stage, set, mean)
            colnames(new) = c("stage", "mean")
            new$sd=aggregate(prop~stage, set, sd)$prop
            new$group=setname
            a_rel=rbind(a_rel, new)


            new = aggregate(number~stage, set, mean)
            colnames(new) = c("stage", "mean")
            new$sd=aggregate(number~stage, set, sd)$number
            new$group=setname
            a_total=rbind(a_total, new)
            return(list(a_total, a_rel))
        }}

        d=read.table("{input}")
        colnames(d) = c("stage", "sample", "sampleName", "number", "lib", "bc")

        d$stage=factor(d$stage, levels=c("raw", "lenFilter", "qualFilter", "winQualFilter", "primerFilter"))

        d$group=NA
        d[d$sample %%in%% c("%(stechlin)s"),]$group="stechlin"
        d[d$sample %%in%% c("%(mock)s"),]$group="mock_community"
        d[d$sample %%in%% c("%(isolate)s"),]$group="isolate"

        at=numeric(0)
        ar=numeric(0)

        env=subset(d, group=="stechlin")
        rv=collectData(env, "stechlin", at, ar)
        at=rv[[1]]
        ar=rv[[2]]


        moc=subset(d, group=="mock_community")
        rv=collectData(moc, "mock community", at, ar)
        at=rv[[1]]
        ar=rv[[2]]

        iso=subset(d, group=="isolate")
        rv=collectData(iso, "isolate samples", at, ar)
        at=rv[[1]]
        ar=rv[[2]]


        #ggplot(subset(ar, stage!="raw"), aes(x=stage, y=mean, fill=group)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + labs(y="proportion of reads remaining")
        ggplot(subset(ar, stage!="raw"), aes(x=stage, y=mean, fill=group)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + labs(y="proportion of reads remaining") + facet_grid(group~.)
        ggsave("{output.graph}", width=16, height=10)

        outtab=dcast(d[!is.na(d$group),], group+sampleName~stage, fun.aggregate=sum, value.var="number")
        outtab=outtab[order(outtab$group, outtab$sampleName),]
        write.table(outtab, "{output.tab}", quote=F, sep="\t", col.names=T, row.names=F)
        """ % {"stechlin": '","'.join([sId for sId, sName in sampleInfo["sampleName"].items() if sName.startswith("STN")]),
               "mock": '","'.join([sId for sId, sType in sampleInfo["sampleType"].items() if sType == "mock"]),
               "isolate": '","'.join([sId for sId, sType in sampleInfo["sampleType"].items() if sType == "isolate"])
              }
        )
        
rule prepPrecluster:
    """Prepare reads for pre-clustering i.e. sort them by mean quality"""
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
    """Pre-cluster at 99% similarity with vsearch, create consensus sequence"""
    input: "preclusters/{sample}_cluInput.fasta"
    output: uc="preclusters/{sample}.uc.txt", cons="consensus/{sample}_consensus.fasta"
    log: "logs/{sample}_pre-cluster.log"
    threads: 6
    shell:
        "%(vsearch)s --usersort --cluster_smallmem {input} --relabel {wildcards.sample}_precluster --sizeout --iddef 0 --id 0.99 --minsl 0.9 --consout {output.cons} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

rule preClusterInfo:
    """Create tables for which read is in which pre cluster and number of reads per pre-cluster"""
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


rule itsx:
    """Run ITSx on pre-cluster"""   
    input: "denovoChimera/{sample}.nochimera.fasta"
    output: "itsx/{sample}.SSU.fasta", "itsx/{sample}.ITS1.fasta", "itsx/{sample}.5_8S.fasta", "itsx/{sample}.ITS2.fasta", "itsx/{sample}.LSU.fasta", "itsx/{sample}.summary.txt", "itsx/{sample}.positions.txt", "itsx/{sample}.full.fasta"
    threads: 6
    log: "logs/{sample}_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 500 -E 1e-4 2> {log}" % config

def sampleMappingInput(wildcards):
    """determine input data for sampleMapping rule according to sampleSet wildcard"""
    if wildcards.sampleSet == "stechlin":
        return expand("itsx/{sampleSet}.full.fasta", sampleSet=stechlin)
    else:
        return expand("itsx/{sampleSet}.full.fasta", sampleSet=samples)

rule getSampleMapping:
    """Create a dictonary (in a pickle file) of which pre-cluster comes from which sample"""
    input: sampleMappingInput
    output: sample="{sampleSet}_preClu2sample.pic"
    run:
        preClu2sample = {}
        for inputFile in input:
            sample = inputFile.rsplit("/", 1)[-1].split(".", 1)[0]
            for rec in SeqIO.parse(open(inputFile), "fasta"):
                preClu2sample[rec.id] = sample
        pickle.dump(preClu2sample, open(output.sample, "wb"))

######### helper functions ##########
def homopoly(kmer):
    """Check if kmer is a homopolymer"""
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
