#!/usr/bin/python2
""" pijpleiding.py  :  pipeliner for the cel-seq protocol

This pipeline accepts one argument only: a configuration file.
Each section in the configuration file is related to a different step in the
pipeline. The options are parsed here, but only few are used directly:

  - pipe_run (boolean): decides whether this step should run at all.
  - pipe_input_files (multiline string): the input file names. Can have several
                                         patterns, each on its own line. Each 
                                         pattern is exapnded to match existing
                                         files.
  - output_dir (string): This directory is created if non existent.

The rest of the parameters are passed as is to the relevant pipe segment.

There are two important constants in this script:
  - SECTIONS - the section names as they appear in the config file
  - SEGMENTS - the functions to be used, in the same order as SECTIONS.

"""


import ConfigParser
import argparse
from glob import glob
import os
from logging import getLogger

import bowtie_wrapper, bc_demultiplex, htseq_wrapper


SECTIONS = ( "bc_demultiplex", "bowtie_wrapper", "htseq_wrapper")
SEGMENTS = ( bc_demultiplex.main, bowtie_wrapper.main, htseq_wrapper.main)

logger = getLogger('pijp')

def main(config_file):
    
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    
    for section, segment in zip(SECTIONS, SEGMENTS):

        if config.getboolean(section, "pipe_run"):
            parameters = dict(config.items(section))
            parameters.pop("pipe_run")
            parameters["input_files"] = []
            for input_glob_pattern in parameters.pop("pipe_input_files").splitlines():
                gl = glob(input_glob_pattern)
                assert (len(gl) != 0 ), "input files not found for pattern "+input_glob_pattern
                parameters["input_files"] += gl

            parameters["input_files"].sort()
            try:
                os.makedirs(parameters['output_dir'])
            except OSError:
                # directory already exists..
                pass
            
            segment(**parameters)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()
    main(args.config_file)
