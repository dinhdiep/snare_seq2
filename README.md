# SNARE-seq2 analysis tools
### Analysis tools for SNARE-seq2 (Scalable Dual-omic Profiling with Single-nucleus Chromatin Accessibility and mRNA Expression Sequencing 2 

### Requirements

Install all dependencies before running. This software have only been tested on a Debian based Linux system but should work on most Linux systems.

The scripts uses Perl, Python3, and R. 

1. samtools (v1.10 or newer): https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
2. pysam: `pip install Cython \ pip install pysam`
3. cutadapt: `python3 -m pip install --user --upgrade cutadapt`
4. deindexer: https://github.com/ws6/deindexer
5. dropEst: https://dropest.readthedocs.io/en/latest/setup
6. STAR: https://github.com/alexdobin/STAR
7. Minimap2: https://github.com/lh3/minimap2 
8. SnapTools: https://github.com/r3fang/snaptool9
9. R packages: httr, broom, igraph, png, rlang, crayon, digest, assertthat, glue, purr, backports, dplyr, tidytext, dropestr, DropletUtils, Seurat, Rcpp, RcppEigen, Rinside, Matrix, optparse, rmarkdown, withr
 
### Download and installation

Download scripts from github:
`git clone http://github.com/dinhdiep/snare_seq2`

There are two main scripts: `SeqToFrag.sh`, `SeqToRDS.sh`. 

RNA count matrices require the separated dT and N6 barcodes to be merged. The R script "create_seurat_kneeplot_2.r" script is provided for this purpose. 

### SeqToFrag.sh
| Input | Explanation of values |
|----------|-------------|
| d | the path to the raw fastq folder |
| f | the fastq id matching the pattern <fastq.id>*R1*fastq.gz, <fastq.id>*R2*fastq.gz, <fastq.id>*R3*fastq.gz |
| l | the library ids for each sample separated by "," (ie sample1,sample2,sample3). Must correspond to the samples in the samples list folder. |
| s | samples list folder that contains the whitelist dT/N6/chromatin barcoded wells for each sample in Round 1 plate |
| a | the path to the minimap2 software folder |
| g | the path to the genome reference fasta file |
| c | the path to the genome reference sequence sizes (ie. <genome>*fai file) |
| n | the genome name used by snaptools (hg38, hg19, mm10, etc). 
| p | path to preseq if sequencing saturation info needed (optional) | 

Outputs: BED formated fragment files, snap files from snaptools, alignment BAM files

Output folders: by_samples_fastq, 01_snaretag, 02_alignment (BAM), 03_snap (count matrices)

### SeqToRDS.sh

| Input | Explanation of values |
|----------|-------------|
| d | the path to the raw fastq folder |
| f | the fastq id matching the pattern <fastq.id>*R1*fastq.gz, <fastq.id>*R2*fastq.gz|
| s | samples list folder that contains the whitelist dT/N6/chromatin barcoded wells for each sample in Round 1 plate |
| p | the path to the dropest software folder |
| a | the path to the STAR software folder |
| i | the path to the STAR index files |
| g | the genes.gtf - the genes annotation file |
| c | the path to cutadapt software folder, (ie ~/.local/bin for python pip installation) |

Outputs: Matrix (.mtx), genes.tsv, barcodes.tsv, R data file with count matrix and dropest diagnostic statistics, alignment BAM files

Output folders: by_samples_fastq, 01_droptag, 02_alignment (BAM), 03_dropEst (count matrices)

Run `create_seurat_kneeplot_2.r 03_dropEst` to generate Seurat rds file with dT and N6 barcodes merged.

### Example usage:

Navigate to the "example_run" directory. Modify `run_atac_rna.sh` shell script to have the correct paths for the variables `softwares_dir` (shared location of deindexer, dropEst, STAR, minimap2, and snaptools) and `ref_dir` (location of all reference folders - see below). 

Go to https://support.10xgenomics.com/single-cell-atac/software/downloads/latest and follow instructions to download GRCh38 (or mm10) reference genome to the `ref_dir` directory.

Go to https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest and follow instructions to download GRCh38 (or mm10) reference genome to the `ref_dir` directory.

Go to https://drive.google.com/drive/folders/1wh6KTd7VqtudZlhEE56JUJ5kd_-_5Irw?usp=sharing to download example fastq files and copy to "example_run/raw_fastq"

To analyse the example fastq files, run: `sh run_atac_rna.sh`

  
