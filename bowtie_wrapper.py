#!/usr/bin/python2
""" run bowtie with specified parameter file 
"""

import subprocess
import logging
import ConfigParser
import os
import shlex
import argparse

logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s %(message)s' , level=logging.INFO)

def build_bowtie_command(fastq_file,  paramfile, output_dir):
    config = ConfigParser.ConfigParser()


    config.read(paramfile)
    index = config.get("bowtie", "index_file")
    threads = config.get("bowtie", "number_of_threads")
    extra_params = config.get("bowtie", "extra_params")

    base_fastq = os.path.splitext(os.path.basename(fastq_file))[0]
    samfile = os.path.join(output_dir, base_fastq + ".sam")
    bowtie_report = os.path.join(output_dir, base_fastq + "_bt2_report.txt")

    bowtie_cmd = "bowtie2 -p {0} {1} {2} -U {3} -S {4} 2> {5}".format(threads, extra_params, index, fastq_file, samfile, bowtie_report)
    return bowtie_cmd


def main(fastq_file, paramfile, output_dir):
    bt2_cmd = build_bowtie_command(fastq_file, paramfile, output_dir)
    
    logging.info("ran  :" + bt2_cmd)
    subprocess.check_call(bt2_cmd, shell=True)
    
    base_fastq = os.path.splitext(os.path.basename(fastq_file))[0]
    file_stem = os.path.join(output_dir, base_fastq )
    sam2bam_cmd = ["samtools","view","-bS","{0}.sam".format(file_stem),"-o","{0}.bam".format(file_stem)]

    logging.info("ran  :" + " ".join(sam2bam_cmd))
    subprocess.check_call(sam2bam_cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastqfile', type=str)
    parser.add_argument('paramfile', type=str)
    parser.add_argument('--output-dir', type=str, default=".")
    args = parser.parse_args()
    main(args.fastqfile, args.paramfile, args.output_dir) 

