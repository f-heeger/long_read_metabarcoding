rule db_getUniteFile:
    output: "%(dbFolder)s/sh_general_release_dynamic_%(uniteVersion)s.fasta" % config
    log: "logs/unite_dl.log"
    shell:
        "cd %(dbFolder)s;"\
        "wget -o {log} -O sh_general_release_dynamic_%(uniteVersion)s.zip %(uniteUrl)s;"\
        "unzip sh_general_release_dynamic_%(uniteVersion)s.zip;"\
        "rm sh_general_release_dynamic_%(uniteVersion)s.zip" % config

rule createUniteTax:
    """create a tsv file with taxon information for UNITE sequences and remove 
    taxon infomration from fasta file"""
    input: "%(dbFolder)s/sh_general_release_dynamic_%(uniteVersion)s.fasta" % config
    output: fasta="%(dbFolder)s/UNITE_%(uniteVersion)s.fasta" % config, tax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/db_createUniteTax.py"

rule creatUniteIndex:
    """Create lambda index file for UNITE data"""
    input: "%(dbFolder)s/UNITE_%(uniteVersion)s.fasta" % config
    output: directory("%(dbFolder)s/UNITE_%(uniteVersion)s.index.lambda" % config)
    threads: 6
    conda:
        "../envs/lambda.yaml"
    shell:
        "lambda2 mkindexn -d {input} -i {output} -t {threads}" % config
        
rule getSilva_main:
    """download SILVA fasta file"""
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz" % config
    shell:
        "cd %(dbFolder)s;" \
        "wget https://www.arb-silva.de/fileadmin/silva_databases/release_%(silvaVersion)s/Exports/SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz" % config

rule getSilva_md5:
    """download silva md5 checksum file"""
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz.md5" % config
    shell: 
        "cd %(dbFolder)s;" \
        "wget https://www.arb-silva.de/fileadmin/silva_databases/release_%(silvaVersion)s/Exports/SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz.md5" % config
    
rule getSilva_test:
    """run md5sum check for silva file"""
    input: gz="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz" % config, md5="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz.md5" % config
    output: touch("%(dbFolder)s/silva_{marker}dl_good" % config)
    shell: 
        "cd %(dbFolder)s;" \
        "md5sum -c SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz.md5" % config

rule unpackSilva:
    """uncompress silva file"""
    input: gz="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta.gz" % config, good="%(dbFolder)s/silva_{marker}dl_good" % config
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta" % config
    shell:
        "cd %(dbFolder)s;" \
        "gunzip SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta.gz; " \
        "touch SILVA_%(silvaVersion)s_{wildcards.marker}Ref_tax_silva_trunc.fasta" % config

rule getSilvaQual:
    """Download silva sequence quality file"""
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality.gz" % config
    shell: "cd %(dbFolder)s;" \
           "wget https://www.arb-silva.de/fileadmin/silva_databases/release_%(silvaVersion)s/Exports/quality/SILVA_%(silvaVersion)s_{wildcards.marker}Ref.quality.gz" % config
           
rule getSilvaQualMd5:
    """Download silva quality md5 checksum file"""
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality.gz.md5" % config
    shell: "cd %(dbFolder)s;" \
           "wget https://www.arb-silva.de/fileadmin/silva_databases/release_%(silvaVersion)s/Exports/quality/SILVA_%(silvaVersion)s_{wildcards.marker}Ref.quality.gz.md5" % config
           
rule getSilvaQual_test:
    """run md5sum check for silva quality file"""
    input: gz="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality.gz" % config, md5="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality.gz.md5" % config
    output: touch("%(dbFolder)s/silva_{marker}dl_qual_good" % config)
    shell: 
        "cd %(dbFolder)s;" \
        "md5sum -c SILVA_%(silvaVersion)s_{wildcards.marker}Ref.quality.gz.md5" % config

rule unpackSilvaQual:
    """uncompress silva quality file"""
    input: gz="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality.gz" % config, good="%(dbFolder)s/silva_{marker}dl_qual_good" % config
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality" % config
    shell:
        "cd %(dbFolder)s;" \
        "gunzip SILVA_%(silvaVersion)s_{wildcards.marker}Ref.quality.gz; " \
        "touch SILVA_%(silvaVersion)s_{wildcards.marker}Ref.quality" % config

rule filterSilva:
    """Filter silva database by sequence quality and chimera probability"""
    input: fasta="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.fasta" % config, qual="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref.quality" % config
    output: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.good.fasta" % config
    log: "%(dbFolder)s/silvaFilter_{marker}.log" % config
    params: minSeqQ=85.0, minPinQ=50.0
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/db_filterSilva.py"

rule getSilvaTaxMap:
    """Download silva taxonomy data"""
    output: "%(dbFolder)s/tax_slv_{marker}_%(silvaVersion)s.txt" % config
    shell:
        "cd %(dbFolder)s;" \
        "wget https://www.arb-silva.de/fileadmin/silva_databases/release_%(silvaVersion)s/Exports/taxonomy/tax_slv_{wildcards.marker}_%(silvaVersion)s.txt" % config

rule getSilvaTaxMapMd5:
    """Download silva taxonomx md5 checksum"""
    output: "%(dbFolder)s/tax_slv_{marker}_%(silvaVersion)s.txt.md5" % config
    shell:
        "cd %(dbFolder)s;" \
        "wget https://www.arb-silva.de/fileadmin/silva_databases/release_%(silvaVersion)s/Exports/taxonomy/tax_slv_{wildcards.marker}_%(silvaVersion)s.txt.md5" % config

rule testSilvaTaxMap:
    """Run silva taxonomy md5sum check"""
    input: md5="%(dbFolder)s/tax_slv_{marker}_%(silvaVersion)s.txt.md5" % config, txt="%(dbFolder)s/tax_slv_ssu_128.txt" % config
    output: touch("%(dbFolder)s/silva_{marker}dl_tax_good" % config)
    shell: 
        "cd %(dbFolder)s;" \
        "md5sum -c tax_slv_{wildcards.marker}_%(silvaVersion)s.txt.md5" % config

def renameSilvaTaxInput(wildcards):
    """determine input data for renameSilvaTax rule (lower case, unlike fasta file)"""
    return "%(dbFolder)s/tax_slv_" % config + wildcards.marker.lower() + "_%(silvaVersion)s.txt" % config #, "%(dbFolder)s/silva_"% config + wildcards.marker.lower() + "dl_tax_good"]

rule renameSilvaTax:
    """Rename silva taxonomy (change marker name to upper case to be consistent 
    with fasta file)"""
    input: renameSilvaTaxInput
    output: "%(dbFolder)s/SILVA_tax_{marker}_%(silvaVersion)s.txt" % config
    shell: "mv {input} {output}"

rule createSlivaTax:
    """Create taxonomy file for silva, only using "canonical" levels"""
    input: fasta="%(dbFolder)                                           s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.good.fasta" % config, tax="%(dbFolder)s/SILVA_tax_{marker}_%(silvaVersion)s.txt" % config
    output: tax="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}_tax.tsv" % config
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/db_createSilvaTax.py"

rule creatSilvaIndex:
    """Create lambda index for silva file"""
    input: "%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.good.fasta" % config
    output: directory("%(dbFolder)s/silva_{marker}_index.lambda" % config)
    threads: 6
    conda:
        "../envs/lambda.yaml"
    shell:
        "lambda2 mkindexn -d {input} -i {output} -t {threads}" % config

rule getRdpLsu:
    """Download rdp LSU file (fasta file and readme with version number)"""
    output: "%(dbFolder)s/current_Fungi_unaligned.fa.gz" % config, "%(dbFolder)s/releaseREADME.txt" % config
    shell: 
        "cd %(dbFolder)s;" \
        "wget http://rdp.cme.msu.edu/download/current_Fungi_unaligned.fa.gz;" \
        "rm releaseREADME.txt;" \
        "wget http://rdp.cme.msu.edu/download/releaseREADME.txt" % config

rule unpackRdpLsu:
    """uncompress rdp data; also check versio number"""
    input: gz="%(dbFolder)s/current_Fungi_unaligned.fa.gz" % config, version="%(dbFolder)s/releaseREADME.txt" % config
    output: "%(dbFolder)s/rdp_LSU_%(rdpVersion)s.fasta" % config
    run:
        v=open(input.version).read().strip()
        if v != config["rdpVersion"]:
            raise RuntimeError("RDP version on server (%s) is different from the one specified in the config file (%s)." % (v, config["rdpVersion"]))
        cmd="cd %(dbFolder)s;" \
        "gunzip current_Fungi_unaligned.fa.gz; " \
        "mv current_Fungi_unaligned.fa rdp_LSU_%(rdpVersion)s.fasta;" \
        "touch rdp_LSU_%(rdpVersion)s.fasta" % config
        shell(cmd)

rule createRdpLsuTax:
    """crete RDP taxonomy file"""
    input: "%(dbFolder)s/rdp_LSU_%(rdpVersion)s.fasta" % config
    output: "%(dbFolder)s/rdp_LSU_%(rdpVersion)s_tax.tsv" % config
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/db_createRdpLsuTax.py"

rule createRdpLsuIndex:
    """create lambda index for RDP database"""
    input: "%(dbFolder)s/rdp_LSU_%(rdpVersion)s.fasta" % config
    output: directory("%(dbFolder)s/rdp_LSU_index.lambda" % config)
    threads: 6
    conda:
        "../envs/lambda.yaml"
    shell:
        "lambda2 mkindexn -d {input} -i {output} -t {threads}"

#rule ssuOverview:
#    """Generate taxonomy level statistics for silva"""
#    input: ssuFasta="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}Ref_tax_silva_trunc.good.fasta" % config, ssuTax="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}_tax.tsv" % config
#    output: ssuOut="%(dbFolder)s/SILVA_%(silvaVersion)s_{marker}_stars.tsv" % config
#    run:
#        ssuSize = {}
#        for rec in SeqIO.parse(open(input.ssuFasta), "fasta"):
#            ssuSize[rec.id] = len(rec)
#        ssuCls = {}
#        for line in open(input.ssuTax):
#            sId, taxStr = line.strip().split("\t")
#            ssuCls[sId] = taxStr.strip(";").split(";")
#        with open(output.ssuOut, "w") as out:
#            for sId, length in ssuSize.items():
#                out.write("%s\t%i\t%s\n" % (sId, length, "\t".join(ssuCls[sId])))

#rule itsOverview:
#    """Generate taxonomy level statistics for UNITE"""
#    input: itsFasta="%(dbFolder)s/UNITE_%(uniteVersion)s.fasta" % config, itsTax="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config
#    output: itsOut="%(dbFolder)s/UNITE_%(uniteVersion)s_stats.tsv" % config
#    run:
#        itsSize = {}
#        for rec in SeqIO.parse(open(input.itsFasta), "fasta"):
#            itsSize[rec.id] = len(rec)
#        itsCls = {}
#        for line in open(input.itsTax):
#            iId, taxStr = line.strip().split("\t")
#            itsCls[iId] = ["Eukaryota"] + [c.split("__", 1)[-1] for c in taxStr.split(";")]
#        with open(output.itsOut, "w") as out:
#            for iId, length in itsSize.items():
#                out.write("%s\t%i\t%s\n" % (iId, length, "\t".join(itsCls[iId])))
#rule lsuOverview:
#    """Generate taxonomy level statistics for RDP LSU"""
#    input: lsuFasta="%(dbFolder)s/rdp_LSU_%(rdpVersion)s.fasta" % config, lsuTax="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_tax.tsv" % config
#    output: lsuOut="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_stats.tsv" % config
#    run:
#        lsuSize = {}
#        for rec in SeqIO.parse(open(input.lsuFasta), "fasta"):
#            lsuSize[rec.id] = len(rec)
#        lsuCls = {}
#        for line in open(input.lsuTax):
#            lId, taxStr = line.strip().split("\t")
#            lsuCls[lId] = taxStr.strip(";").split(";")
#        with open(output.lsuOut, "w") as out:
#            for lId, length in lsuSize.items():
#                out.write("%s\t%i\t%s\n" % (lId, length, "\t".join(lsuCls[lId])))

#rule dbOverviewPlot:
#    """plot taxonomy level statistics for all thress databases"""
#    input: ssu="%(dbFolder)s/SILVA_%(silvaVersion)s_SSU_stars.tsv" % config, its="%(dbFolder)s/UNITE_%(uniteVersion)s_stats.tsv" % config, lsu="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_stats.tsv" % config
#    output: length="dbOverview_len.pdf", tax="dbOverview_tax.pdf"
#    conda:
#        "envs/ggplot.yaml"
#    script:
#        "scripts/plotDbOverview.R"
#        
#rule taxComp:
#    input: ssu="%(dbFolder)s/SILVA_%(silvaVersion)s_SSU_tax.tsv" % config, its="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config, lsu="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_tax.tsv" % config
#    output: "dbMissingGroups.tsv"
#    run:
#        ssuCls = {}
#        for line in open(input.ssu):
#            sId, taxStr = line.strip().split("\t")
#            t_ssuCls = taxStr.strip(";").split(";")
#            cur = ssuCls
#            for entry in t_ssuCls:
#                try:
#                    cur = cur[entry]
#                except KeyError:
#                    cur[entry] = {}
#                    cur = cur[entry]
#        itsCls = {}
#        for line in open(input.its):
#            iId, taxStr = line.strip().split("\t")
#            t_itsCls = ["Eukaryota"] + [c.split("__", 1)[-1] for c in taxStr.split(";")]
#            cur = itsCls
#            for entry in t_itsCls:
#                try:
#                    cur = cur[entry]
#                except KeyError:
#                    cur[entry] = {}
#                    cur = cur[entry]
#        lsuCls = {}
#        for line in open(input.lsu):
#            lId, taxStr = line.strip().split("\t")
#            t_lsuCls = taxStr.strip(";").split(";")
#            cur = lsuCls
#            for entry in t_lsuCls:
#                try:
#                    cur = cur[entry]
#                except KeyError:
#                    cur[entry] = {}
#                    cur = cur[entry]
#        tax = {"ssu": ssuCls, "its": itsCls, "lsu": lsuCls}
#        miss = {}
#        #SSU
#        for line in open(input.ssu):
#            sId, taxStr = line.strip().split("\t")
#            t_ssuCls = taxStr.strip(";").split(";")
#            if t_ssuCls[1] != "Fungi":
#                continue
#            for mrk1, mrk2 in [("ssu","lsu"), ("ssu", "its")]:
#                cur = tax[mrk2]
#                for r, sEntry in enumerate(t_ssuCls):
#                    try:
#                        cur = cur[sEntry]
#                    except KeyError:
#                        break
#                try:
#                    miss[(mrk1, mrk2, sEntry, str(r), t_ssuCls[2])].append(sId)
#                except KeyError:
#                    miss[(mrk1, mrk2, sEntry, str(r), t_ssuCls[2])] = [sId]
#        #ITS
#        for line in open(input.its):
#            iId, taxStr = line.strip().split("\t")
#            t_itsCls = ["Eukaryota"] + [c.split("__", 1)[-1] for c in taxStr.split(";")]
#            if t_itsCls[1] != "Fungi":
#                continue
#            for mrk1, mrk2 in [("its", "ssu"), ("its", "lsu")]:
#                cur = tax[mrk2]
#                for r, iEntry in enumerate(t_itsCls):
#                    try:
#                        cur = cur[iEntry]
#                    except KeyError:
#                        break
#                try:
#                    miss[(mrk1, mrk2, iEntry, str(r), t_itsCls[2])].append(iId)
#                except KeyError:
#                    miss[(mrk1, mrk2, iEntry, str(r), t_itsCls[2])] = [iId]

#        #LSU
#        for line in open(input.lsu):
#            lId, taxStr = line.strip().split("\t")
#            t_lsuCls = taxStr.strip(";").split(";")
#            if t_lsuCls[1] != "Fungi":
#                continue
#            for mrk1, mrk2 in [("lsu", "ssu"), ("lsu", "its")]:
#                cur = tax[mrk2]
#                for r, lEntry in enumerate(t_lsuCls):
#                    try:
#                        cur = cur[lEntry]
#                    except KeyError:
#                        break
#                try:
#                    miss[(mrk1, mrk2, lEntry, str(r), t_lsuCls[2])].append(lId)
#                except KeyError:
#                    miss[(mrk1, mrk2, lEntry, str(r), t_lsuCls[2])] = [lId]
#        with open(output[0], "w") as out:
#            for key, value in miss.items():
#                out.write("%s\t%i\n" % ("\t".join(key), len(value)))
#                
#rule mockDbOccurence:
#    input: ssu="%(dbFolder)s/SILVA_%(silvaVersion)s_SSU_tax.tsv" % config, its="%(dbFolder)s/UNITE_%(uniteVersion)s_tax.tsv" % config, lsu="%(dbFolder)s/rdp_LSU_%(rdpVersion)s_tax.tsv" % config
#    output: "mockOccurence.tsv"
#    run:
#        mockList=["Clavariopsis aquatica",
#              "Clonostachys rosea",
#              "Trichoderma reesei",
#              "Ustilago maydis",
#              "Saccharomyces cerevisiae",
#              "Mortierella elongata",
#              "Cystobasidium laryngis",
#              "Metschnikowia reukaufii",
#              "Penicillium piscarium",
#              "Exobasidium vaccinii",
#              "Phanerochaete chrysosporium",
#              "Leucosporidium scottii",
#              "Penicillium brevicompactum",
#              "Davidiella tassiana"]
#        with open(output[0], "w") as out:
#            for mock in mockList:
#                for line in shell("set +e; grep -F \"%s\" %s; echo 0 > /dev/null" % (mock.replace(" ", "_"), input.ssu), iterable=True):
#                    out.write("SSU\t%s\t%s\n" % (mock, line.strip()))
#                for line in shell("set +e; grep -F \"%s\" %s; echo 0 > /dev/null" % (mock.replace(" ", "_"), input.its), iterable=True):
#                    out.write("ITS\t%s\t%s\n" % (mock, line.strip()))
#                for line in shell("set +e; grep -F \"%s\" %s; echo 0 > /dev/null" % (mock, input.lsu), iterable=True):
#                    out.write("LSU\t%s\t%s\n" % (mock, line.strip()))
