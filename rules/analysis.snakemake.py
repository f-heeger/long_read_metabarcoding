import pickle

rule concatItsxResult:
    """Concatenate ITSx results from different samples"""
    input: expand("itsx/{sample}.{{marker}}.fasta", sample=samples)
    output: "catItsx/all.{marker}.fasta"
    shell:
        "cat {input} > {output}"

def catIsolatesInput(wildcards):
    """determine input data for catIsolate rule"""
    return ["itsx/%s.%s.fasta" % (s, wildcards.marker) for s in isolates[wildcards.spec]]

rule concatIsolates:
    """Concatenate ITSx results for an isolate sample"""
    input: catIsolatesInput
    output:"isolates/{spec}.{marker}.fasta"
    shell:
        "cat {input} > {output}"

def otuInput(wildcards):
    """determine input file for OTU clustering according to sample wildcard"""
    if wildcards.sample == "all":
        return "catItsx/all.full.fasta"
#    elif wildcards.sample in isolates:
#        return "isolates/%s.full.fasta" % wildcards.sample
    else:
        return "itsx/%s.full.fasta" % wildcards.sample

rule otuCluster:
    """Cluster OTUs by ITS sequence with vsearch"""
    input: otuInput
    output: fasta="otus/{sample}_{ident}otus.fasta", uc="otus/{sample}_{ident}otus.uc.tsv"
    log: "logs/{sample}_{ident}otuClustering.log"
    threads: 3
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --cluster_size {input} --relabel otu --sizein --sizeout --iddef 0 --id 0.{wildcards.ident} --minsl 0.9 --centroids {output.fasta} --uc {output.uc} --threads {threads} --log {log} &> /dev/null" % config

def otuReadsInput(wildcards):
    """determine input data for otuReads rule according to sample wildcard"""
    if wildcards.sample == "all":
        return ["otus/all_%sotus.uc.tsv" % wildcards.ident] + expand("preclusters/{sample}_cluInfo.tsv", sample=samples)
#    elif wildcards.sample in isolates:
#        return ["otus/%s_%sotus.uc.tsv" % (wildcards.sample, wildcards.ident)] + expand("preclusters/{sample}_cluInfo.tsv", sample=isolates[wildcards.sample])
    else:
        return  ["otus/%s_%sotus.uc.tsv" % (wildcards.sample, wildcards.ident), "preclusters/%s_cluInfo.tsv" % (wildcards.sample)]

rule otuReads:
    """Generate different infos for the OTUs (as dictonaries in pickle files)
    
    *otus.size.tsv: number of reads in the OTU
    *otuReadInfo.pic: reads for each OTU
    *otu_preClusterInfo.pic: preClusters for each OTU
    *otu_repSeq.pic: representative sequence of the OTU
    """
    input: otuReadsInput
    output: size="otus/{sample}_{ident}otus.size.tsv", info="otus/{sample}_{ident}otuReadInfo.pic", preInfo="otus/{sample}_{ident}otu_preClusterInfo.pic", repInfo="otus/{sample}_{ident}otu_repSeq.pic"
    script:
        "../scripts/analysis_otuReads.py"


def transferOtusInput(wildcards):
    """determine input data for transferOtus rule according to sample wildcard"""
    if wildcards.sample == "all":
        return ["otus/all_%sotu_repSeq.pic" % wildcards.ident, "catItsx/all.%s.fasta" % wildcards.marker]
#    elif wildcards.sample in isolates:
#        return ["otus/%s_%sotu_repSeq.pic" % (wildcards.sample, wildcards.ident), "isolates/%s.%s.fasta" % (wildcards.sample, wildcards.marker)]
    else:
        return ["otus/%s_%sotu_repSeq.pic" % (wildcards.sample, wildcards.ident), "itsx/%s.%s.fasta" % (wildcards.sample, wildcards.marker)]


rule transferOtus:
    """Create file with non-ITS marker sequence (SSU or LSU) per OTU"""
    input: transferOtusInput
    output: "otus/{sample}_{ident}otus_{marker}.fasta"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/analysis_transferOtus.py"

rule alignToUnite:
    """Search for local sequence matches of ITS sequences in the UNITE database with lambda"""
    input: clu="otus/{sample}_{ident}otus.fasta", db="%(dbFolder)s/UNITE_%(uniteVersion)s.index.lambda" % config
    output: "lambda/{sample}.{ident}otu_vs_UNITE.m8"
    log: "logs/{sample}_{ident}otu_lambda.log"
    threads: 6
    conda:
        "../envs/lambda.yaml"
    shell:
        "lambda2 searchn -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -t {threads} &> {log}" % config

rule classifyITS:
    """Classify OTU based on ITS matches with LCA approach"""
    input: lam="lambda/{sample}.{ident}otu_vs_UNITE.m8", tax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
    output: "taxonomy/{sample}_{ident}otu_ITS.class.tsv"
    log: "logs/{sample}_{ident}otuClass.log", "logs/{sample}_{ident}otu_itsTax.log"
    script:
        "../scripts/analysis_classifyITS.py"

rule alignToSilva:
    """Search for local sequence matches of SSU sequences in the SILVA database with lambda"""
    input: clu="otus/{sample}_{ident}otus_SSU.fasta", db="%(dbFolder)s/silva_SSU_index.lambda" % config
    output: "lambda/{sample}.{ident}otu_SSU_vs_SILVA.m8"
    log: "logs/{sample}_{ident}otu_SSU_lambda.log"
    threads: 6
    conda:
        "../envs/lambda.yaml"
    shell:
        "lambda2 searchn -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -n 10000 -t {threads} -b -2 -x 30 --adaptive-seeding F &> {log}" % config

rule classifySSU:
    """Classify OTU based on SSU matches with LCA approach"""
    input: lam="lambda/{sample}.{ident}otu_SSU_vs_SILVA.m8", tax="%(dbFolder)s/SILVA_%(silvaVersion)s_SSU_tax.tsv" % config
    output: "taxonomy/{sample}_{ident}otu_SSU.class.tsv"
    log: "logs/{sample}_SSU_{ident}otuClass.log", "logs/{sample}_{ident}otu_SSU_tax.log", "logs/{sample}_{ident}otu_SSU_tiling.log"
    conda:
        "../envs/networkx.yaml"
    script:
        "../scripts/analysis_classifySSU.py"

rule alignToRdp:
    """Search for local sequence matches of LSU sequences in the RDP LSU database with lambda"""
    input: clu="otus/{sample}_{ident}otus_LSU.fasta", db="%(dbFolder)s/rdp_LSU_index.lambda" % config
    output: "lambda/{sample}.{ident}otu_LSU_vs_RDP.m8"
    log: "logs/{sample}_{ident}otuLSU_lambda.log"
    threads: 6
    conda:
        "../envs/lambda.yaml"
    shell:
        "lambda2 searchn -q {input.clu} -i {input.db} -o {output} --output-columns \"std qlen slen\" -n 5000 -t {threads} -b -2 -x 40 --adaptive-seeding F &> {log}" % config

rule classifyLSU:
    """Classify OTU based on LSU matches with LCA approach"""
    input: lam="lambda/{sample}.{ident}otu_LSU_vs_RDP.m8", tax="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_tax.tsv" % config
    output: "taxonomy/{sample}_{ident}otu_LSU.class.tsv"
    log: "logs/{sample}_LSU_{ident}otuClass.log", "logs/{sample}_{ident}otu_LSU_tax.log", "logs/{sample}_{ident}otu_LSU_tiling.log"
    conda:
        "../envs/networkx.yaml"
    script:
        "../scripts/analysis_classifyLSU.py"

rule combineCls:
    """Create table with OTU classifications with SSU, ITS, LSU"""
    input: ssu="taxonomy/{sample}_{ident}otu_SSU.class.tsv", its="taxonomy/{sample}_{ident}otu_ITS.class.tsv", lsu="taxonomy/{sample}_{ident}otu_LSU.class.tsv", size="otus/{sample}_{ident}otus.size.tsv"
    output: "taxonomy/{sample}_{ident}_comb.class.tsv"
    run:
        ssu = {}
        for line in open(input.ssu):
            ssuId, ssuTax = line.strip().split("\t")
            ssu[ssuId] = ssuTax
        lsu = {}
        for line in open(input.lsu):
            lsuId, lsuTax = line.strip().split("\t")
            lsu[lsuId] = lsuTax
        its = {}
        for line in open(input.its):
            itsName, itsTax = line.strip().split("\t")
            itsId, itsSizeStr = itsName.strip(";").split(";")
            its[itsId] = itsTax
        with open(output[0], "w") as out:
            for line in open(input.size):
                itsId, sizeStr = line.strip().split("\t")
                size = int(sizeStr)
                itsTax = its.get(itsId, "unknown")
                lsuTax = lsu.get(itsId, "unknown")
                ssuTax = ssu.get(itsId, "unknown")
                out.write("%s\t%i\t%s\t%s\t%s\n" % (itsId, size, ssuTax, itsTax, lsuTax))

rule createOtuTab:
    """Create a OTU table, with size, classiifcation and occurence information"""
    input: tax="taxonomy/{sampleSet}_{ident}_comb.class.tsv", otu2pClu="otus/{sampleSet}_{ident}otu_preClusterInfo.pic", sample="{sampleSet}_preClu2sample.pic"
    output: "{sampleSet}_otu{ident}_table.tsv"
    run:
        otu2pClu = pickle.load(open(input.otu2pClu, "rb"))
        pCluSample = pickle.load(open(input.sample, "rb"))
        sampleOrder = list(set(pCluSample.values()))
        tab = [["otu", "totalSize", "ssuTax", "itsTax", "lsuTax"] + sampleOrder]
        for line in open(input.tax):
            oId, size, ssuTax, itsTax, lsuTax = line.strip().split("\t")
            pClusters = {}
            for pClu in otu2pClu[oId]:
                tSample = pCluSample[pClu]
                try:
                    pClusters[tSample].append(pClu)
                except KeyError:
                    pClusters[tSample] = [pClu]
            tab.append([oId, size, ssuTax, itsTax, lsuTax])
            for sample in sampleOrder:
                inThisSample = 0
                for pClu in pClusters.get(sample, []):
                    inThisSample += int(pClu.split(";")[1].split("=")[1])
                tab[-1].append(str(inThisSample))
        with open(output[0], "w") as out:
            for line in tab:
                out.write("\t".join(line)+"\n")

rule collectClsStats:
    """Create table of classifications (including depth) as input for ggplot"""
    input: cls="taxonomy/{sampleSet}_97_comb.class.tsv"
    output: complete="taxonomy/{sampleSet}_97_comb.stats.tsv" 
    run:
        ranks = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        tab = []
        for line in open(input.cls):
            oId, size, ssuCls, itsCls, lsuCls = line.strip().split("\t")
            tab.append((oId, int(size), ssuCls, itsCls, lsuCls))

        tab.sort(key=lambda x: x[1], reverse=True)

        with open(output.complete, "w") as out:
            for oId, size, ssuClsStr, itsClsStr, lsuClsStr in tab:
                ssuDepth = 0
                ssuCls = dict(zip(ranks, ["NA"]*8))
                if ssuClsStr != "unknown":
                    ssuArr = ssuClsStr.strip(";").split(";")
                    ssuDepth = len(ssuArr)-1
                    for r, rank in enumerate(ranks):
                        try:
                            ssuCls[rank] = ssuArr[r]
                        except IndexError:
                            break
                lsuDepth = 0
                lsuCls = dict(zip(ranks, ["NA"]*8))
                if lsuClsStr != "unknown":
                    lsuArr = lsuClsStr.strip(";").split(";")
                    lsuDepth = len(lsuArr)-1
                    for r, rank in enumerate(ranks):
                        try:
                            lsuCls[rank] = lsuArr[r]
                        except IndexError:
                            break
                itsDepth = 0
                itsCls = dict(zip(ranks, ["NA"]*8))
                if itsClsStr != "unknown":
                    itsArr = itsClsStr.strip(";").split(";")
                    itsDepth = len(itsArr)
                    itsCls["kingdom"] = "Eukaryota"
                    for r, rank in enumerate(ranks[1:]):
                        try:
                            itsCls[rank] = itsArr[r][3:]
                        except IndexError:
                            break
                out.write("%s\t%i\tssu\t%s\t%i\n" % (oId, size, "\t".join([ssuCls[r] for r in ranks]), ssuDepth))
                out.write("%s\t%i\tits\t%s\t%i\n" % (oId, size, "\t".join([itsCls[r] for r in ranks]), itsDepth))
                out.write("%s\t%i\tlsu\t%s\t%i\n" % (oId, size, "\t".join([lsuCls[r] for r in ranks]), lsuDepth))

rule clsSummary:
    """Collect some summary stats of how many OTUs were assigend to taxonomic 
    ranks for the paper abstract"""
    input: "taxonomy/{sampleSet}_97_comb.stats.tsv"
    output: "taxonomy/{sampleSet}_97_clsStats.tsv"
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/clsSummary.R"

rule plotClsComp:
    """Create plots of classifications depth"""
    input: all="taxonomy/{sampleSet}_97_comb.stats.tsv"
    output: depth="{sampleSet}_clsComp_depth.pdf", depthFungi="{sampleSet}_clsComp_depth_fungi.pdf", block="{sampleSet}_clsComp_basic.pdf"
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/plotClsSummary.R"

rule compareCls:
    """Compare classifications from different markers for plotting with ggplot
    
    Is a OTU that is classified by marker A at a certain level also classified 
    by marker B and C at this level.
    """
    input: cls="taxonomy/{sampleSet}_97_comb.class.tsv"
    output: diff="{sampleSet}_clsDiff.tsv", comp="{sampleSet}_clsComp.tsv", diffStat="{sampleSet}_clsDiffStat.tsv"
    script:
        "../scripts/analysis_compareClass.py"

rule plotDiff:
    """Plot comparison data of classification with different markers"""
    input: "{sampleSet}_clsDiffStat.tsv"
    output: bars="{sampleSet}_clsDiffStat.svg", prop="{sampleSet}_clsDiffPropSame.svg"
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/plotDiff.R"

#def plotChimeraEnvInput(wildcards):
#    """determine input data for plotChimera rule according to sample wildcard"""
#    if wildcards.sample == "all":
#        return ["denovoChimera/%s.chimeraReport.tsv" % s for s in samples]
#    elif wildcards.sample in isolates:
#        return ["denovoChimera/%s.chimeraReport.tsv" % s for s in isoaltes[wildcards.sample]]
#    else:
#        return "denovoChimera/%s.chimeraReport.tsv" % wildcards.sample

#rule plotChimeraEnv:
#    """Plot chimera data for environmental samples"""
#    input: plotChimeraEnvInput
#    output: tab="denovoChimera/{sample}_chimeraTable.tsv", relative="denovoChimera/{sample}_chimerasRelativeBarplot.svg"
#    conda:
#        "../envs/ggplot.yaml"
#    script:
#        "../scripts/plotChimeraEnv.py"



