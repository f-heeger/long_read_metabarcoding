export PATH=/srv/public/shared/software/miniconda3/bin:$PATH


source activate snakemake5.4.2

snakemake --use-conda $@

conda deactivate
