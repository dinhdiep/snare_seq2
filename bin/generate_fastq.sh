#!/bin/bash

usage() { echo "Usage: $0 [-d <run.directory>] [-m <string>] [-s SampleSheet.csv] [-c <num.cores>] [-b <barcodes.mismatch.num>] [-l <lanes: ie 1-4 ] [-o <out.dir>]" 1>&2; exit 1; }

while getopts ":d:m:s:c:b:l:o:" options; do
    case $options in
    	d ) run_dir=$OPTARG;;
        m ) bases_mask_string=$OPTARG;;
	s ) sample_sheet=$OPTARG;;
	c ) num_cores=$OPTARG;;
	b ) mismatch_num=$OPTARG;;
	l ) lanes=$OPTARG;;
	o ) output_dir=$OPTARG;;

    esac
done
shift $(($OPTIND - 1))

if [ -z "${run_dir}" ] || [ -z "${bases_mask_string}" ] || [ -z "${sample_sheet}" ]; then
    usage
fi

run_info_xml="$run_dir/RunInfo.xml"


### Setting up to generate FASTQ ###
cp $run_info_xml $output_dir/

### Generate Undetermined Fastq File ###
bcl2fastq -R $run_dir --output $output_dir --create-fastq-for-index-reads --use-bases-mask $bases_mask_string --sample-sheet $sample_sheet --mask-short-adapter-reads=8 --minimum-trimmed-read-length=8 --barcode-mismatches $mismatch_num --no-lane-splitting -p $num_cores --tiles s_[$lanes]

