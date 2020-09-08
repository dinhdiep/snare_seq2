import pysam as ps
import sys
from collections import Counter
from itertools import product
import subprocess as sp
import multiprocessing as mp
import os
import time
from functools import partial


def main():
    read1_fastq_file = sys.argv[1]
    read2_fastq_file = sys.argv[2]
    out_prefix = sys.argv[3]
    index_cfg = sys.argv[4]
    index_list = sys.argv[5]

    ## Time the script
    start = time.time()

    ## Hardcoded parameters
    #index_list = "config/R1_barcode_list_split-seq"
    n_mismatch = 2

    ## Split out index files and generate config 
    print "Splitting index file..."
    split_index(read2_fastq_file, out_prefix + ".Round1", 86, 94)
     
    ## Run deindexer
    print "Running fastq deindexing..."
    run_deindexer(out_prefix, read1_fastq_file, read2_fastq_file, index_list, index_cfg, n_mismatch)

    print "Time elapsed:\t" + str(time.time() - start)
    



def run_deindexer(out_file_prefix, read1_fastq_file, read2_fastq_file, index_list, index_cfg, n_mismatch = 1):
    if not os.path.isdir(out_file_prefix + "_deindexed_fastq"):
        os.mkdir(out_file_prefix + "_deindexed_fastq")
    
    deindexer_cmd = "deindexer -f rrb -c " + index_cfg + " -b " + index_list + \
        " -o " + out_file_prefix + "_deindexed_fastq" + \
        " -m " + str(n_mismatch) + \
        " " + read1_fastq_file + " " + \
        " " + read2_fastq_file + " " + \
        out_file_prefix + ".Round1.idx.fastq "

    print deindexer_cmd 
    sp.call("ulimit -n 500000; " + deindexer_cmd, shell = True)

 	


def split_index(idx_file, out_file_prefix, start_idx, stop_idx):
    f_in = ps.FastxFile(idx_file)
    f_out = open(out_file_prefix + ".idx.fastq", mode = "w")
    for read in f_in:
        read.sequence = read.sequence[start_idx:stop_idx]
        read.quality = read.quality[start_idx:stop_idx]
        f_out.write(str(read) + "\n")
    f_in.close()
    f_out.close()


if __name__ == "__main__": main()
