softwares_dir=$HOME
ref_dir=$HOME

../bin/SeqToRDS.sh -d raw_fastq -f N713_S1 -s plate1_lists \
   -p $softwares_dir/dropEst \
   -a $softwares_dir/STAR-2.5.1b/bin/Linux_x86_64 \
   -i $ref_dir/refdata-cellranger-GRCh38-3.0.0/star \
   -g $ref_dir/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
   -c .local/bin/ 

../bin/SeqToFrag.sh -d raw_fastq -f N521_S1 -l sample1 -s plate1_lists \
   -a $HOME/minimap2/ \
   -g $HOME/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa \
   -c $HOME/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa.fai \
   -n hg38
