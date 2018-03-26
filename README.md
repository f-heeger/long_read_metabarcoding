# Long Read Metabarcoding
This pipeline was developed to analyze PacBio CSS amplicon data of the fungi rRNA operon, but might also be useful for other eukaryotes and with different primers.

The results of this analysis are described in our manuscript: [*Long-read DNA metabarcoding of ribosomal rRNA in the analysis of fungi from aquatic environments*](https://www.biorxiv.org/content/early/2018/03/15/283127).

The exact version used for the results in teh manuscript can be found as release v1.0.

## Dependecies
The pipeline is implemented as [snakemake](http://snakemake.readthedocs.io/en/stable/) work flow. It contains a mix of rules directly written in python and rules that call external tools. The path for external tools can be defined in the config.json file.

The workflow uses several python packages and external tools that have to be installed. They were tested with the version given below, but newer versions might work as well. In addition three different database are used as reference, that will be downloaded automatically in the version given in the config file.

### Python pakages
* [snakemake](http://snakemake.readthedocs.io/en/stable/) (version 3.5.5)
* [biopython](https://pypi.python.org/pypi/biopython/1.70) (version 1.66)
* [networkx](https://pypi.python.org/pypi/networkx/) (version 1.11)

### External Tools
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.2)
* [multiqc](http://multiqc.info/) (version 0.6)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/) (version 1.9.1)
* [vsearch](https://github.com/torognes/vsearch) (version 2.4.3)
* [lambda](https://seqan.github.io/lambda/) (version 1.9.2)
* [itsx](http://microbiology.se/software/itsx/) (version 1.0.11)
* [blasr](https://github.com/PacificBiosciences/blasr)   (github version from [Oct 13, 2016](https://github.com/PacificBiosciences/blasr/tree/16b158dd09d79fc6274da11fa154c45f0f689ef0))
* [SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) (version 2.9.0)
### Databases
* [UNITE](https://unite.ut.ee/) (20.11.2016)
* [SILVA](https://www.arb-silva.de/) (128)
* [RDP LSU](http://rdp.cme.msu.edu/misc/resources.jsp) (11.5)


## Running the analysis
The pipeline is setup to reproduce the analysis in the manuscript, but can be adapted to work on other data as well.

After installing depnedencies and configuring paths in the config file (`config.json`), the analysis can be run with `snakemake -s run_analysis.snakemake.py` (the other snakemake files are included in this one). You might also want to use `-j` to give multiple processors to snakemake. 

The following resulting files were used in the manusript (with some adjustment for readability):

* Figure 2: `chimeraCyclesRelativeBarplot.svg`
* Figure 3: was created from `mock/clusterGraph/Lib4-0018_clusterGraphCls.tsv`, `mock/clusterGraph/Lib4-0018_clusterGraphEdges.tsv` and `mock/clusterGraph/Lib4-0018_clusterGraphClsLab.tsv` using Cytoscape
* Figure 4: `all_clsComp_depth_fungi.svg`
* Figure 5: `all_clsDiffStat.svg`

Interesting graphs about error rates and read assignment that were not used in the manuscript can be found in the `mapping` folder.


## Rules

The following is a list of the main rules in the workflow. Many simple rules are not specifically documented here.

### getFullCls
For each pre-cluster representative sequences get a classification. This can be either  i) CHIMERA if the reference base chimera detection called this as chimeric (Y) or possibly chimeric (?),  ii) UNKNOWN if the sequence was not called as chimeric, but not match to an isolate consensus sequence was found or iii) the name of species, this is given by the highest scoring match to a isolate consensus sequence for non-chimeric sequences.

### fullMapping
Run blasr to map representatives of non-chimeric pre-clusters against the isolate consensus sequences. The following parameters are used: -m 5 to get tabular output, --bestn 50 to get a maximum of 50 hits for each query and --minPctSimilarity 90 to get only hits with at least 90% identity.

### removeChimeraRef
Run vsearch to remove chimeras with reference based approach. The --uchime_ref parameter is used to run the reference based chimera detection algorithm and the --db parameter to give the isolate consensus sequences selected by the getFullRef rule as reference sequences.

### getFullRef
For each isolate sample get the consensus sequence of the biggest pre-cluster. Will give an error if there are more than one pre-cluster with 10 or more reads for one sample. Compares sequences for replicates of each species (if available) and writes a warning to the log file if a difference is encountered.

### getCorrectCls
Get “correct” classifications for reads in each OTU according to mappings to isolate consensus sequences. Each OTU might have multiple species with read numbers listed.

### classifyLSU
Classify OTUs by matches of LSU sequences to the RDP LSU database. See main methods section for details and parameters.

### alignToRdp
Run lambda to get local alignments of each OTU LSU representative to the RDP LSU database. Lambda is run with the following parameters: --output-columns "std qlen slen" to get query length and subject length along with the default columns in the output table (these are used for coverage computation later), -p blastn to run in blastn (nucleotide vs nucleotide) mode, -nm 5000 to get more matches per query sequence (this is important due to the high number of similar sequences in the database), -b -2 to th square root of the query length as width for the banded alignment optimization (because higher indel rate and, even more important, uneven insertion/deletion ratio cause the alignment to leave the band), -x 40 to reduce the x-drop value which can cause alignments to be terminated prematurely and -as F to disable adaptive seeding which normally is used to reduce number of hits. Parameters were optimized to allow for alignments for all mock community species.

### classifySSU
Classify OTUs by matches of SSU sequences to the SILVA database. See main methods section for details and parameters.

### alignToSilva
Run lambda to get local alignments of each OTU SSU representative to the SILVA database. Lambda is run with the following parameters: --output-columns "std qlen slen" to get query length and subject length along with the default columns in the output table (these are used for coverage computation later), -p blastn to run in blastn (nucleotide vs nucleotide) mode, -nm 20000 to get more matches per query sequence (this is important due to the high number of similar sequences in the database), -b -2 to th square root of the query length as width for the banded alignment optimization (because higher indel rate and, even more important, uneven insertion/deletion ratio cause the alignment to leave the band), -x 30 to reduce the x-drop value which can cause alignments to be terminated prematurely and -as F to disable adaptive seeding which normally is used to reduce number of hits. Parameters were optimized to allow for alignments for all mock community species.

### classifyITS
Classify OTUs by matches to the UNITE database. See main methods section for details and parameters.

### alignToUnite
Run lambda to get local alignments of each OTU representative to the UNITE database. Lambda is run with default parameters except for: --output-columns "std qlen slen" to get query length and subject length along with the default columns in the output table (these are used for coverage computation later) and -p blastn to run in blastn (nucleotide vs nucleotide) mode. 

### otuCluster
Run vsearch to cluster OTUs at 97% identity threshold. The following parameters are used for vsearch: --cluster_size to choose cluster seeds by descending pre-cluster size (according to size annotation), --relabel otu to name OTUs with out and running number instead of using the first sequence as a name, --sizein --sizeout to read and write size annotation, --iddef 0 to use the identity definition of CD-Hit (see vsearch manual), --id 0.97 to use 97% identity threshold and --minsl 0.9 to only accept alignments of at least 90% coverage for similarity computation. In addition --centroids is used to output a representative centroid sequence for each OTU.

### itsx
Run ITSx to separate the different regions of the rRNA operon. The following parameters were set for ITSx: -t . to use HMM models from all available taxonomic groups, --save_regions SSU,ITS1,5.8S,ITS2,LSU to save separate files for all different regions of the rRNA operon, --complement F to not allow reverse-complement detection (all sequences were orientated in forward direction in the primerFilter rule), --partial 500 to allow for partial matches (we do not have complete SSU and LSU sequences in the amplicon) and -E 1e-4 to allow for HMM hits with slightly lower e-values (this was optimized to make sure that all rRNA operons in the mock community species were recognized).

### removeChimera
Run vsearch to remove chimeras in de novo mode. Vsearch is run with the –uchime_denovo parameter. All parameters for the chimera detection algorithm are left at default values. Vsearch automatically uses size annotations form the pre-clustering step for its greedy algorithm.

### preCluster
Run vsearch to create pre-clusters at 99% identity threshold. The following parameters are used for vsearch: --usersort --cluster_smallmem to choose cluster seeds in the order the sequences are sorted in the input file, --relabel {wildcards.sample}_precluster to name pre-clusters according to the given scheme instead of using the first sequence as a name, --sizeout to add size annotation to the output, --iddef 0 to use the identity definition of CD-Hit (see vsearch manual), --id 0.99 to use 99% identity threshold and --minsl 0.9 to only accept alignments of at least 90% coverage for similarity computation. In addition --consout is used to generate consensus sequences for each pre-cluster.

### prepPrecluster
Reads are sorted by descending mean quality (see qualityFilter rule for computation). This helps to use high quality reads as cluster seeds for pre-clusters in the next step.

### filterPrimer
Primers are found and cut with cutadapt. Cutadapt is run with default parameters except for --trimmed-only to only retain reads were the primer was found and -O 10.  Cutadapt is configured to search for both the forward and the reverse primer at the start of the sequence. For sequences where the forward primer was found, the reverse-complemented reverse primer is search at the end of the sequence with an additional run of cutadapt. Accordingly for sequences where the reverse primer was found, the reverse-complemented forward primer is searched at the end of the sequence with another run of cutadapt. In the end sequences with forward-reverse primer combination and reverse-complemented sequences with reverse-forward primer combination are concatenated into one file. 

### windowQualFilter
For overlapping windows of size 8 the mean error rate is computed from the Phred scores with the same formula as in the qualityFilter rule (except that S is the substring in the windiw instead of the whole sequence). If any window in a sequence has a mean error rate of 0.9 or higher the sequence is removed. 

### qualityFilter
Mean error rate per sequence is computed from the Phred score given in the fastq file. Sequences with an error rate of 0.4% or more are written to a separate file (not further used).

### lengthFilter
Sequences with length above 6,500 or below 3,000 are printed to separate files (not used further).

### filterSilva
Filter SILVA sequences by the quality and pintail (chimera probability) values given in the database. Only sequences with a quality value of at least 85 and pintail value of at least 50 are retained.
