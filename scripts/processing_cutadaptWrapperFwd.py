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
    comp = {"A": "T", "T": "A", "C": "G", "G": "C",
            "W": "W", "M": "K", "K": "M", "R": "Y",
            "Y": "R", "S": "S", "H": "D", "D": "H",
            "V": "B", "B": "V", "N": "N"}
    
    fBC = snakemake.config["samples"][snakemake.wildcards.sample]["fwdBarcodeSeq"]
    rBC = snakemake.config["samples"][snakemake.wildcards.sample]["revBarcodeSeq"]
    rc_rBC = "".join(comp[base] for base in rBC[::-1])
    
    cmd = ["cutadapt", "-g" , "%s...%s" % (fBC, rc_rBC), "--discard-untrimmed", "--cores", str(snakemake.threads), "--error-rate", str(snakemake.config["primerErr"]), "-o", snakemake.output[0], snakemake.input[0]]
    print(cmd)
    exit(subprocess.call(cmd, stdout=open(snakemake.log[0], "w")))
