import re

rule createBam:
    input: "raw/{part}.bax.h5"
    output: "raw/{part}.subreads.bam"
    conda:
        "../envs/bax2bam.yaml"
    shell:
        "bax2bam --subread -o raw/{wildcards.part} {input}"

rule ccs:
    input: "raw/{part}.subreads.bam"
    output: "ccs/{part}.ccs.bam"
    log: "logs/ccs/{part}.ccs.log", "logs/ccs/{part}_css_report.txt"
    threads: 6
    conda:
        "../envs/ccs.yaml"
    shell:
        "ccs --minPasses %(minPass)s --minPredictedAccuracy %(minAcc)s -j {threads} --logFile {log[0]} --reportFile {log[1]} {input} {output}" % config


rule mergeCcs:
    input: expand("ccs/{part}.ccs.bam", part=config["parts"])
    output: "ccs/all.ccs.bam"
    threads: 6
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge -@{threads} {output} {input}"

rule barcodeFiles:
    output: "barcodes.fasta"
    script:
        "../scripts/processing_createBarcodeFile.py"


rule demultiplex:
    input: bam="ccs/all.ccs.bam", bc="barcodes.fasta"
    output: expand("fastq/{sample}.bam", sample=samples.keys())
    conda:
        "../envs/lima.yaml"
    script:
        "../scripts/processing_demultiplex.py"


def createFastqGzInput(wildcards):
    sraId = config["samples"][wildcards.sample]["sraId"]
    if len(sraId) == 0:
        return ["fastq/%s.bam" % wildcards.sample]
    elif re.match("SRR[0-9]{7}", sraId) is None:
        raise RuntimeError("Data source %s for sample %s is invalid" % ())
    else:
        return []

rule createFastqGz:
    input: createFastqGzInput
    output: "fastq/{sample}.fastq.gz"
    conda:
        "../envs/createFastqGz.yaml"
    script:
        "../scripts/createFastqGz.py"

rule fastqc:
    """Run fastqc an all samples"""
    input: "fastq/{sample}.fastq.gz"
    output: "QC/{sample}_fastqc.html"
    threads: 6
    conda:
        "../envs/fastqc.yaml"
    shell: 
        "cp {input} QC/{wildcards.sample}.fastq; fastqc -o QC -t {threads} QC/{wildcards.sample}.fastq; rm QC/{wildcards.sample}.fastq" % config

rule multiqc:
    """Combine QC reports into one with multiqc"""
    input: expand("QC/{sample}_fastqc.html", sample=samples)
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc -f --interactive -o QC QC/*_fastqc.zip" % config

rule concatFastq:
    input: expand("fastq/{sample}.fastq.gz", sample=samples)
    output: "fastq/all.fastq.gz"
    shell:
        "zcat {input} | gzip -c > {output}"

rule lengthFilter:
    """Filter reads by minimum and maximum length"""
    input: "fastq/{sample}.fastq.gz"
    output: right="lenFilter/{sample}_rightLen.fastq.gz", long="lenFilter/{sample}_tooLong.fastq.gz", short="lenFilter/{sample}_tooShort.fastq.gz"
    log: "logs/{sample}_lenFilter.log"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_lengthFilter.py"

rule qualityFilter:
    """Filter reads by minimum mean quality"""
    input: fastq="lenFilter/{sample}_rightLen.fastq.gz"
    output: good="qualFilter/{sample}_goodQual.fastq.gz", bad="qualFilter/{sample}_badQual.fastq.gz", info="qualFilter/{sample}_qualInfo.tsv"
    log: "logs/{sample}_qualityFilter.log"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_qualityFilter.py"

#rule qualityVsLength:
#    """Generate data for plotting length vs quality"""
#    input: fasta="lenFilter/{sample}_rightLen.fastq", qual="qualFilter/{sample}_qualInfo.tsv"
#    output: tsv="qualFilter/{sample}_LenVsQual.tsv"
#    run:
#        qual = {}
#        for line in open(input.qual):
#            rId, q = line.strip().split("\t")
#            qual[rId] = float(q)
#        with open(output.tsv, "w") as out:
#            for rec in SeqIO.parse(open(input.fasta), "fastq"):
#                out.write("%s\t%s\t%i\t%f\n" % (wildcards.sample, rec.id, len(rec), qual[rec.id]))

#rule concatQualVsLength:
#    """Concatenate quality vs length data for all samples"""
#    input: expand("qualFilter/{sample}_LenVsQual.tsv", sample=samples)
#    output: "qualFilter/lenVsQual.tsv"
#    shell:
#        "cat {input} > {output}"

#rule plotQualVsLength:
#    """Plot quality vs length"""
#    input: "qualFilter/lenVsQual.tsv"
#    output: "lenVsQual.png"
#    conda:
#        "envs/ggplot.yaml"
#    script:
#        "scripts/qualVsLength.R"
        

rule windowQualFilter:
    """Filter by sliding window mean quality"""
    input: fastq="qualFilter/{sample}_goodQual.fastq.gz"
    output: good="windowQualFilter/{sample}_goodQual.fasta", stat="windowQualFilter/{sample}_stat.tsv"
    log: "logs/{sample}_winQualityFilter.log"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_windowQualityFilter.py"

rule catWindowQual:
    """concatenate window quality data for all environmental samples"""
    input: expand("windowQualFilter/{sample}_stat.tsv", sample=samples)
    output: "windowQualFilter/envStats.tsv"
    shell:
        "cat {input} > {output}"


rule filterPrimer53:
    """Filter sequences by occurence of primer sequences and cut primer sequences
    with forward primer in the front and reverse primer in the end"""
    input: "windowQualFilter/{sample}_goodQual.fasta"
    output: fastq=temp("primers/{sample}_53.fasta")
    log: "logs/{sample}_53_primer.log"
    threads: 6
    conda:
        "../envs/cutadapt.yaml"
    script:
        "../scripts/processing_cutadaptWrapperFwd.py"

rule filterPrimer35:
    """Filter sequences by occurence of primer sequences and cut primer sequences
    with reverse complement reverse-primer in the front and 
    reverse complement forward-primer in the end"""
    input: "windowQualFilter/{sample}_goodQual.fasta"
    output: fastq=temp("primers/{sample}_35.fasta")
    log: "logs/{sample}_35_primer.log"
    threads: 6
    conda:
        "../envs/cutadapt.yaml"
    script:
        "../scripts/processing_cutadaptWrapperRev.py"
        

rule combinFilterPrimer:
    """combine the forward and reverse Primer filtered sequences"""
    input: fwd="primers/{sample}_53.fasta", rev="primers/{sample}_35.fasta"
    output: fastq="primers/{sample}_primer.fasta"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_combineFilterPrimer.py"
    


def readNumbers_rawInput(wildcards):
    return [s["path"] for s in config["samples"].values()]

rule readNumbers_raw:
    """Count raw reads"""
    input: expand("fastq/{sample}.fastq.gz", sample=samples)
    output: "readNumbers/rawReadNumbers.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_readNumbersRaw.py"


rule readNumbers_maxLen:
    """Count reads after length filter"""
    input: length=expand("lenFilter/{sample}_rightLen.fastq.gz", sample=samples)
    output: "readNumbers/lenReadNumbers.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_readNumbersMaxLen.py"

rule readNumbers_qual:
    """Count reads after average quality filter"""
    input: qual=expand("qualFilter/{sample}_goodQual.fastq.gz", sample=samples)
    output: "readNumbers/qualReadNumbers.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_readNumbersQual.py"

rule readNumbers_winQual:
    """Countr reads after sliding window quality filter"""
    input: qual=expand("windowQualFilter/{sample}_goodQual.fasta", sample=samples)
    output: "readNumbers/winQualReadNumbers.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_readNumbersWinQual.py"

rule readNumbers_primer:
    """Count reads after primer filtering"""
    input: primer=expand("primers/{sample}_primer.fasta", sample=samples)
    output: "readNumbers/primerReadNumbers.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_readNumbersPrimer.py"

rule catReadNumber:
    """Concatenate all read numbers for plotting"""
    input: "readNumbers/rawReadNumbers.tsv", "readNumbers/lenReadNumbers.tsv", "readNumbers/qualReadNumbers.tsv", "readNumbers/winQualReadNumbers.tsv", "readNumbers/primerReadNumbers.tsv"
    output: "readNumbers/readNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            for inFile in input:
                for line in open(inFile):
                    stage, sample, number = line.strip().split("\t")
                    out.write("%s\t%s\t%s\t%s\n" % (stage, sample, samples[sample]["name"], number))

rule plotReadNumber:
    """Plot read number after each filtering step"""
    input: "readNumbers/readNumbers.tsv", "samples.tsv"
    output: filtering="readNumbers/readNumbersFiltering.svg", groups="readNumbers/readNumbersGroups.svg"
    script:
       "../scripts/plotReadNumber.R"
        
rule prepPrecluster:
    """Prepare reads for pre-clustering i.e. sort them by mean quality"""
    input: fasta="primers/{sample}_primer.fasta", qual="qualFilter/{sample}_qualInfo.tsv"
    output: "preclusters/{sample}_cluInput.fasta"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_prepPreCluster.py"

rule preCluster:
    """Pre-cluster at 99% similarity with vsearch, create consensus sequence"""
    input: "preclusters/{sample}_cluInput.fasta"
    output: uc="preclusters/{sample}.uc.txt", cons="consensus/{sample}_consensus.fasta"
    log: "logs/{sample}_pre-cluster.log"
    threads: 6
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --usersort --cluster_smallmem {input} --relabel {wildcards.sample}_precluster --sizeout --iddef 0 --id 0.99 --minsl 0.9 --consout {output.cons} --uc {output.uc} --threads {threads} --log {log} > /dev/null" % config

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

rule removeChimeraDenovo:
    input: seqs="consensus/{sample}_consensus.fasta"
    output: fasta="denovoChimera/{sample}.nochimera.fasta", tsv="denovoChimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_denovoChimera.log"
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --uchime_denovo {input.seqs} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config


rule itsx:
    """Run ITSx on pre-cluster"""   
    input: "denovoChimera/{sample}.nochimera.fasta"
    output: "itsx/{sample}.SSU.fasta", "itsx/{sample}.ITS1.fasta", "itsx/{sample}.5_8S.fasta", "itsx/{sample}.ITS2.fasta", "itsx/{sample}.LSU.fasta", "itsx/{sample}.summary.txt", "itsx/{sample}.positions.txt", "itsx/{sample}.full.fasta"
    threads: 6
    log: "logs/{sample}_itsx.log"
    conda:
        "../envs/itsx.yaml"
    shell:
        "ITSx -t . -i {input} -o itsx/{wildcards.sample} --save_regions SSU,ITS1,5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 100 -E 1e-4 2> {log}" % config

rule collectLogData:
    input: raw="readNumbers/rawReadNumbers.tsv", length=expand("logs/{sample}_lenFilter.log", sample=samples), qual=expand("logs/{sample}_qualityFilter.log", sample=samples), window=expand("logs/{sample}_winQualityFilter.log", sample=samples), primer1=expand("logs/{sample}_53_primer.log", sample=samples), primer2=expand("logs/{sample}_35_primer.log", sample=samples), chimera=expand("logs/{sample}_denovoChimera.log", sample=samples), itsx=expand("itsx/{sample}.summary.txt", sample=samples)
    output: "logs/logData.tsv"
    script:
        "../scripts/processing_collectLogData.py"

rule plotRemoval:
    input: "logs/logData.tsv"
    output: "readsAfterFiltering.pdf"
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/processing_plotRemoval.R"

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
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/processing_getSampleMapping.py"

