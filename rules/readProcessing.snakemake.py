
#rule getData:
#    """Download read data from the SRA"""
#    output: "%(inFolder)s/{sample}.fastq" % config
#    run:
#        sId = sampleInfo["sraID"][wildcards.sample]
#        shell("%(fastq-dump)s" % config + " %s -Z > {output}" % sId)

rule fastqc:
    """Run fastqc an all samples"""
    input: "%(inFolder)s/{sample}.fastq" % config
    output: "QC/{sample}_fastqc.html"
    threads: 6
    conda:
        "envs/fastqc.yaml"
    shell: 
        "fastqc -o QC -t {threads} {input}" % config

rule multiqc:
    """Combine QC reports into one with multiqc"""
    input: expand("QC/{sample}_fastqc.html", sample=allSamples)
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -f --interactive -o QC QC/*_fastqc.zip" % config

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
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/qualVsLength.R"
        

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
    threads: 6
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -a ^%(forward_primer)s...%(reverse_primer)s$ --discard-untrimmed --cores {threads} --error-rate %(primerErr)f -o {output} {input} > {log}" % config


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

#rule plotReadNumber:
#    """Plot read number after each filtering step"""
#    input: "readNumbers/readNumbers.tsv"
#    output: graph="readNumbers.svg", tab="readNumbersTab.tsv"
#    script:
#       "scripts/readProcessing.R"
        
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
    conda:
        "envs/vsearch.yaml"
    shell:
        "vsearch --usersort --cluster_smallmem {input} --relabel {wildcards.sample}_precluster --sizeout --iddef 0 --id 0.99 --minsl 0.9 --consout {output.cons} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

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
    conda:
        "envs/itsx.ymal"
    shell:
        "ITSx -t . -i {input} -o itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 100 -E 1e-4 2> {log}" % config

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
