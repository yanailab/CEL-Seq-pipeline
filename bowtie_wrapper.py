#!/usr/bin/python2
""" run bowtie with specified parameter file 
"""

import subprocess
from logging import getLogger
import os
import shlex
import argparse
from glob import glob

logger = getLogger('pijp.bowtie_wrapper')

def build_bowtie_command(fastq_file,  index_file, number_of_threads, output_dir, extra_params):


    base_fastq = os.path.splitext(os.path.basename(fastq_file))[0]
    samfile = os.path.join(output_dir, base_fastq + ".sam")
    bowtie_report = os.path.join(output_dir, base_fastq + "_bt2_report.txt")

    bowtie_cmd = "bowtie2 -p {0} {1} {2} -U {3} -S {4} 2> {5}".format(number_of_threads, extra_params, index_file, fastq_file, samfile, bowtie_report)
    return bowtie_cmd


def main(input_files, index_file, number_of_threads, output_dir, extra_params):


    for fastq_file in input_files:
        bt2_cmd = build_bowtie_command(fastq_file, index_file, number_of_threads, output_dir, extra_params)
    
        logger.info("ran : " + bt2_cmd)
        subprocess.check_call(bt2_cmd, shell=True)
    
