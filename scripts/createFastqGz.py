import subprocess

sraId = snakemake.config["samples"][snakemake.wildcards.sample]["sraId"]

if len(sraId) == 0:
    subprocess.run(["bam2fastq", "-o", "fastq/%s" % snakemake.wildcards.sample, snakemake.input[0]])
else:
    subprocess.run(["fastq-dump", sraId, "-Z", "--gzip"], stdout=open(snakemake.output[0], "w"))
