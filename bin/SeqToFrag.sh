#!/bin/bash

usage() { echo "Usage: $0 [-d <raw.fastq.directory>] [-f <fastq.id>] [-l <library.id>] [-s <samples.list.directory>] [-a <minimap2.path>] [-g <genome.fa>] [-c <genome.fai>] [-n <genome.name> ] [ -p path to preseq ]" 1>&2; exit 1; }

while getopts ":d:f:s:a:g:c:n:p:b:l:" options; do
  case $options in
    d ) raw_fastq_dir=$OPTARG;;
    f ) fastq_id=$OPTARG;;
    s ) list_dir=$OPTARG;;
    a ) minimap2_path=$OPTARG;;
    g ) genome_fasta=$OPTARG;;
    c ) genome_fai=$OPTARG;;
    n ) genome_name=$OPTARG;;
    p ) preseq_path=$OPTARG;;
    b ) blacklist_path=$OPTARG;;
    l ) library_id=$OPTARG;;
  esac
done
shift $(($OPTIND - 1))

if [ -z "${raw_fastq_dir}" ] || [ -z "${fastq_id}" ] || [ -z "${list_dir}" ] || [ -z "${minimap2_path}" ] || [ -z "${genome_fasta}" ] || [ -z "${genome_name}" ] || [ -z "${genome_fai}" ] ; then
  usage
fi

snare2_bin=`dirname "$(realpath $0)"`
snare2_config="$snare2_bin/../config"
cur_dir=`pwd`

##### Demultiplex with deindexer #####
raw_R1=`find $raw_fastq_dir -name '*_R1*fastq.gz' -print | grep $fastq_id`
raw_R2=`find $raw_fastq_dir -name '*_R2*fastq.gz' -print | grep $fastq_id`
raw_R3=`find $raw_fastq_dir -name '*_R3*fastq.gz' -print | grep $fastq_id`
python $snare2_bin/demultiplex_SNARE2.py \
  $raw_R1 $raw_R2 $raw_R3 $fastq_id $snare2_config/R1_barcode_cfg $snare2_config/R1_barcode_list_SNARE2


sample_files=`ls $list_dir/* | cat | sort`

##### Merge barcode 1 #####

mkdir by_samples_fastq

for list in $sample_files
do
  sample_list=`basename "$list"`
  echo $sample_list
  sample=`echo $sample_list | sed 's/.list//g'`
  rm -f by_samples_fastq/${sample}.${fastq_id}_R1.fastq
  rm -f by_samples_fastq/${sample}.${fastq_id}_R2.fastq
  rm -f by_samples_fastq/${sample}.${fastq_id}_R3.fastq
  rm -f by_samples_fastq/${sample}.${fastq_id}_R1.fastq.gz
  rm -f by_samples_fastq/${sample}.${fastq_id}_R2.fastq.gz
  rm -f by_samples_fastq/${sample}.${fastq_id}_R3.fastq.gz
  for barcode in `cat $list`
  do
    cat ${fastq_id}_deindexed_fastq/*_${barcode}_R1.fastq >> by_samples_fastq/${sample}.${fastq_id}_R1.fastq
    cat ${fastq_id}_deindexed_fastq/*_${barcode}_R2.fastq >> by_samples_fastq/${sample}.${fastq_id}_R2.fastq
    cat ${fastq_id}_deindexed_fastq/*_${barcode}_R3.fastq >> by_samples_fastq/${sample}.${fastq_id}_R3.fastq
  done
done

rm -fr ${fastq_id}_deindexed_fastq


##### SNARETag #####
cd $cur_dir
mkdir 01_snareTag
cd $cur_dir/01_snareTag
i=0
for list in $sample_files
do
  sample_list=`basename "$list"`
  sample=`echo $sample_list | sed 's/.list//g'`
  file_id="${sample}.${fastq_id}"
  if [ -z "${library_id}" ]
  then
    cur_id="${file_id}"
  else
    ids_list=($(echo $library_id | tr "," "\n"))
    cur_id=${ids_list[$i]}
    echo "library id given"
  fi
  $snare2_bin/get_barcodes_by_pos_SNARE2.pl $cur_dir/by_samples_fastq/${file_id}_R1.fastq $cur_dir/by_samples_fastq/${file_id}_R2.fastq $cur_dir/by_samples_fastq/${file_id}_R3.fastq ${cur_id}
  mv ${cur_id}.R1.fastq ${file_id}.R1.fastq
  mv ${cur_id}.R3.fastq ${file_id}.R3.fastq
  pigz -p 12 $cur_dir/by_samples_fastq/${file_id}_R1.fastq
  pigz -p 12 $cur_dir/by_samples_fastq/${file_id}_R2.fastq
  pigz -p 12 $cur_dir/by_samples_fastq/${file_id}_R3.fastq
  ((i++))
done

##### snapTools align-paired-end #####
cd $cur_dir
mkdir 02_alignment
cd $cur_dir/02_alignment

for list in $sample_files
do
  sample_list=`basename "$list"`
  sample=`echo $sample_list | sed 's/.list//g'`
  file_id="${sample}.${fastq_id}"
  snaptools align-paired-end \
    --input-reference=$genome_fasta \
    --input-fastq1=../01_snareTag/${file_id}.R1.fastq \
    --input-fastq2=../01_snareTag/${file_id}.R3.fastq \
    --output-bam=${file_id}.bam \
    --aligner=minimap2 \
    --path-to-aligner $minimap2_path \
    --read-fastq-command=cat --min-cov 0 --num-threads 12 --if-sort True --tmp-folder ./ --overwrite True
  rm -f ../01_snareTag/${file_id}.R1.fastq
  rm -f ../01_snareTag/${file_id}.R3.fastq
  if [ -z "${preseq_path}" ]
  then 
    echo "no preseq"
  else
    $preseq_path/bam2mr ${file_id}.bam | sort -k1,1d -k2,2n > ${file_id}.bam.mr
    sed 's/:/\t/g' ${file_id}.bam.mr | sed 's/FRAG\t//g' | awk '{print $1"_"$2"_"$4"_"$5}' | $snare2_bin/tabulateList.pl | cut -f 2 | $snare2_bin/tabulateList.pl | sort > ${file_id}.cov
    $preseq_path/preseq lc_extrap -o ${file_id}.yield.txt -H ${file_id}.cov
  fi
  if [ -z "${blacklist_path}" ]
  then
    echo "no blacklist"
  else
    bedtools intersect -v -abam ${file_id}.bam -b $blacklist_path > ${file_id}.rmsk.bam
  fi
done

##### snapTools snap-pre #####
cd $cur_dir
mkdir 03_snap
cd $cur_dir/03_snap

for list in $sample_files
do
  sample_list=`basename "$list"`
  sample=`echo $sample_list | sed 's/.list//g'`
  file_id="${sample}.${fastq_id}"
  if [ -z "${blacklist_path}" ]
  then
    echo "no blacklist"
  else
    file_id="${sample}.${fastq_id}.rmsk" 
  fi
  $snare2_bin/collapseBarcodesBam_2.pl $snare2_config/whitelist ../02_alignment/${file_id}.bam ../02_alignment/${file_id}.collapsed 12
  $snare2_bin/getFragments.pl ../02_alignment/${file_id}.collapsed.bam > ${file_id}.nsorted.bed 
  snaptools snap-pre  \
    --input-file=${file_id}.nsorted.bed \
    --output-snap=${file_id}.collapsed.snap  \
    --genome-name=$genome_name  \
    --genome-size=$genome_fai \
    --min-mapq=30  \
    --min-flen=0  \
    --max-flen=1000  \
    --keep-chrm=TRUE  \
    --keep-single=TRUE  \
    --keep-secondary=FALSE  \
    --overwrite=True  \
    --min-cov=100 \
    --verbose=FALSE
  snaptools snap-add-bmat --snap-file ${file_id}.collapsed.snap --bin-size-list 1000 5000 10000 --verbose TRUE
done

