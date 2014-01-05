#!/usr/bin/python2
""" run bowtie with specified parameter file 
"""

import subprocess
from logging import getLogger
import os
import shlex

logger = getLogger('pijp.scythe_wrapper')


def main(input_files, output_dir, adaptor_file, extra_params=""):


    for fastq_file in input_files:
        base_fastq = os.path.basename(fastq_file)
        output_file =  os.path.join(output_dir, base_fastq)
        # sometimes base_fastq contains gz, while the output is plaintext
        if output_file.endswith('.gz'):
            output_file = output_file[:-3]

        scythe_cmd = "scythe -q sanger {3} -a {0} {1} -o {2}".format(adaptor_file, fastq_file, output_file, extra_params)
        logger.info("ran  :" + scythe_cmd)
        pro = subprocess.Popen(scythe_cmd, shell=True, stderr = subprocess.PIPE)
        (x, stderr) = pro.communicate()
        assert (pro.returncode == 0 ), "scythe error %d : %s" % (pro.returncode, stderr)
    
        logger.info("Printing scythe output:")
        logger.info(stderr)
