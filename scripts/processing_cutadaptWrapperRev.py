import subprocess
import os

#this script checks if the input file is empty. In this case cutadapt would 
# throw and error (see: https://github.com/marcelm/cutadapt/issues/407) so 
# it just ctreates an empty file as output. If the input file has content 
# cutadapt is called.

if os.path.getsize(snakemake.input[0]) == 0:
    open(snakemake.output[0], "w").close()
    exit(0)
else:
    cmd = ["cutadapt", "-g" , "%(rev_primer)s...%(rv_fwd_primer)s" % snakemake.config, "--discard-untrimmed", "--cores",  snakemake.threads, "--error-rate", snakemake.config["primerErr"], "-o", snakemake.output[0], snakemake.input[0]]

    exit(subprocess.call(cmd, stdout=open(snakemake.log[0], "w")))
