
rule removeChimeraRef:
    input: seqs="primers/{sample}_primer.fasta", ref="isolateSeqs.fasta"
    output: fasta="refChimera/{sample}.nochimera.fasta", tsv="refChimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_refChimera.log"
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --uchime_ref {input.seqs} --db {input.ref} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule subset:
    input: seqs="primers/{sample}_primer.fasta", tsv="refChimera/{sample}.chimeraReport.tsv"
    output: "manualCheck/{sample}_subset.fasta"
    conda:
        "../envs/biopython.yaml"
    params: N=100
    script:
        "../scripts/chimera_subset.py"

rule sampleNamesTab:
    output: "mockSamples.tsv"
    run:
        with open(output[0], "w") as out:
            for sId, sample in sampleInfo["sampleName"].items():
                if sampleInfo["sampleType"][sId] == "mock":
                    out.write("%s\t%s\n" % (sId, sample))

#rule plotChimera:
#    input: "mockSamples.tsv", expand("refChimera/{sample}.chimeraReport.tsv", sample=mockSamples)
#    output: total="chimeraTotalBarplot.svg", tab="chimeraTable.tsv", relative_cy="chimeraCyclesRelativeBarplot.svg", relative_in="chimeraInputRelativeBarplot.svg", subset_cy="chimeraCyclesSubset.svg", subset_in="chimeraInputSubset.svg"
#    params: subset=200
#    conda:
#        "envs/ggplot.yaml"
#    script:
#        "scripts/plotChimera.R"

rule removeChimeraDenovo:
    input: seqs="consensus/{sample}_consensus.fasta"
    output: fasta="denovoChimera/{sample}.nochimera.fasta", tsv="denovoChimera/{sample}.chimeraReport.tsv"
    log: "logs/{sample}_denovoChimera.log"
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --uchime_denovo {input.seqs} --nonchimeras {output.fasta} --uchimeout {output.tsv} &> {log}" % config

rule compareChimera:
    input: denovo="denovoChimera/{sample}.chimeraReport.tsv", ref="refChimera/{sample}.chimeraReport.tsv", cluster2read="preclusters/{sample}_cluInfo.tsv"
    output: "refChimera_{sample}_comp.tsv"
    run:
        clu2read = {}
        for line in open(input.cluster2read):
            clu, read = line.strip().split("\t")
            try:
                clu2read[clu].append(read)
            except KeyError:
                clu2read[clu] = [read]
        ref={}
        for line in open(input.ref):
            arr = line.strip().split("\t")
            scr=arr[0]
            readId = arr[1]
            cls = arr[-1]
            ref[readId] = (cls, scr)
        with open(output[0], "w") as out:
            out.write("readId\tclusterNr\tclusterId\tdenovoScr\tdenovo\trefScr\tref\n")
            for l,line in enumerate(open(input.denovo)):
                arr = line.strip().split("\t")
                scr=arr[0]
                cluId = arr[1].split("=", 1)[1].split(";", 1)[0]
                cls = arr[-1]
                for rId in clu2read[cluId]:
                    out.write("%s\tcluster%i\t%s\t%s\t%s\t%s\t%s\n" % (rId, l, cluId, scr, cls, ref[rId][1], ref[rId][0]))


rule plotChimeraComp:
    input: "refChimera_Lib4-0018_comp.tsv"
    output: "chimera_comp_sankey.svg"
    script:
        "../scripts/plotChimeraComp.R"
