#!/bin/bash

usage() { echo "Usage: $0 [-d <raw.fastq.directory>] [-f <fastq.id>] [-s <samples.list.directory>] [-p <dropest.path>] [-a <star.path>] [-i <star.index>] [-g <gtf.file> ] [-c cutadapt.path]" 1>&2; exit 1; }

while getopts ":d:f:s:p:a:i:g:c:" options; do
  case $options in
    d ) raw_fastq_dir=$OPTARG;;
    f ) fastq_id=$OPTARG;;
    s ) list_dir=$OPTARG;;
    p ) dropest_path=$OPTARG;;
    a ) star_path=$OPTARG;;
    i ) star_index=$OPTARG;;
    g ) gtf_file=$OPTARG;;
    c ) cutadapt_dir=$OPTARG;;
  esac
done
shift $(($OPTIND - 1))

if [ -z "${raw_fastq_dir}" ] || [ -z "${fastq_id}" ] || [ -z "${list_dir}" ] || [ -z "${dropest_path}" ] || [ -z "${star_path}" ] || [ -z "${star_index}" ] || [ -z "${gtf_file}" ]  || [ -z "${cutadapt_dir}" ] ; then
  usage
fi

dropest_dir="$dropest_path/build"
config_file="$dropest_path/config/split_seq.xml"
snare2_bin=`dirname "$(realpath $0)"`
snare2_config="$snare2_bin/../config"
droptag_config_file="$snare2_config/split_seq.xml"
cur_dir=`pwd`

##### Demultiplex with deindexer #####
raw_R1=`find $raw_fastq_dir -name '*_R1*fastq.gz' -print | grep $fastq_id`
raw_R2=`find $raw_fastq_dir -name '*_R2*fastq.gz' -print | grep $fastq_id`

python $snare2_bin/demultiplex_split-seq.py \
	$raw_R1 $raw_R2 $fastq_id $snare2_config/R1_barcode_cfg $snare2_config/R1_barcode_list_split-seq

sample_files=`ls $list_dir/* | cat`

##### Merge barcode 1 #####

mkdir by_samples_fastq
for list in $sample_files
do
  sample_list=`basename "$list"`
  echo $sample_list
  sample=`echo $sample_list | sed 's/.list//g'`
  rm -f by_samples_fastq/${sample}.${fastq_id}_R1.fastq
  rm -f by_samples_fastq/${sample}.${fastq_id}_R2.fastq
  rm -f by_samples_fastq/${sample}.${fastq_id}_R1.fastq.gz
  rm -f by_samples_fastq/${sample}.${fastq_id}_R2.fastq.gz
  for barcode in `cat $list`
  do
    cat ${fastq_id}_deindexed_fastq/*_${barcode}_R1.fastq >> by_samples_fastq/${sample}.${fastq_id}_R1.fastq
    cat ${fastq_id}_deindexed_fastq/*_${barcode}_R2.fastq >> by_samples_fastq/${sample}.${fastq_id}_R2.fastq
  done
done

rm -fr ${fastq_id}_deindexed_fastq

##### DropTag #####
cd $cur_dir
mkdir 01_dropTag
cd $cur_dir/01_dropTag

for list in $sample_files
do
  sample_list=`basename "$list"`
  sample=`echo $sample_list | sed 's/.list//g'`
  file_id="${sample}.${fastq_id}"
  $cutadapt_dir/cutadapt --no-indels --cores=12 --overlap 8 -e 0.2 -a ACGTACTGCAX -o ${file_id}.trimmed.R2.fastq -p ${file_id}.trimmed.R1.fastq --discard-trimmed \
     $cur_dir/by_samples_fastq/${file_id}_R2.fastq $cur_dir/by_samples_fastq/${file_id}_R1.fastq
  $dropest_dir/droptag -p 6 -c $droptag_config_file -l ${file_id} ${file_id}.trimmed.R2.fastq ${file_id}.trimmed.R1.fastq
  pigz -p 4 $cur_dir/by_samples_fastq/${file_id}_R1.fastq
  pigz -p 4 $cur_dir/by_samples_fastq/${file_id}_R2.fastq
  rm ${file_id}.trimmed.R2.fastq
  rm ${file_id}.trimmed.R1.fastq
done

##### STAR #####
cd $cur_dir
mkdir 02_alignment
cd $cur_dir/02_alignment

for list in $sample_files
do
  sample_list=`basename "$list"`
  sample=`echo $sample_list | sed 's/.list//g'`
  file_id="${sample}.${fastq_id}"
  $star_path/STAR  --runThreadN 12 --limitOutSJcollapsed 3000000 --genomeDir $star_index --readFilesCommand zcat --outSAMtype BAM Unsorted --readFilesIn ../01_dropTag/${file_id}.trimmed.R1.fastq.tagged.1.fastq.gz --outFileNamePrefix ${file_id}.
  mv ${file_id}.Aligned.out.bam ${file_id}.out.bam
done

##### DropEst #####
cd $cur_dir
mkdir 03_dropEst
cd $cur_dir/03_dropEst
for list in $sample_files
do
  sample_list=`basename "$list"`
  sample=`echo $sample_list | sed 's/.list//g'`
  file_id="${sample}.${fastq_id}"
  $dropest_dir/dropest -w -M -u -G 20 -L iIeEBA -u -m -F -g $gtf_file -c $droptag_config_file -o ${file_id} -l ${file_id} $cur_dir/02_alignment/${file_id}.out.bam 
done
