import os

from snakemake.shell import shell

shell("lima --different --ccs --split-bam-named %s %s fastq/sample.bam" % (snakemake.input.bam, snakemake.input.bc))

for sId, sample in snakemake.config["samples"].items():
    baseSrcFwd = "fastq/sample.%(fwdBarcodeId)s--%(revBarcodeId)s" % sample
    baseSrcRev = "fastq/sample.%(revBarcodeId)s--%(fwdBarcodeId)s" % sample
    baseTrg = "fastq/%s" % sId
    for ext in ["bam", "bam.pbi", "subreadset.xml"]:
        if os.path.exists("%s.%s" % (baseSrcFwd, ext)):
            os.replace("%s.%s" % (baseSrcFwd, ext), "%s.%s" % (baseTrg, ext))
        elif os.path.exists("%s.%s" % (baseSrcRev, ext)):
            os.replace("%s.%s" % (baseSrcRev, ext), "%s.%s" % (baseTrg, ext))
        else:
            print("No file for a combination of %s-%s "
                  "with extension %s exists" % (sample["fwdBarcodeId"], 
                    sample["revBarcodeId"], ext))
            open("%s.%s" % (baseTrg, ext), "w").close()
        
