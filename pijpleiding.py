#!/usr/bin/python2
""" pijpleiding.py  :  pipeliner for yanailab 

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
import logging
import json

import bowtie_wrapper, bc_demultiplex, htseq_wrapper, clean_up


SECTIONS = ( "bc_demultiplex", "bowtie_wrapper", "htseq_wrapper", "clean_up")
SEGMENTS = ( bc_demultiplex.main, bowtie_wrapper.main, htseq_wrapper.main, clean_up.main)

LOGFORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logger = logging.getLogger('pijp')
logging.basicConfig(level=logging.INFO, format = LOGFORMAT)
log_formatter = logging.Formatter(LOGFORMAT)

def main(config_file):
    
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    logger.info("===== Started pijpleiding with config file : %s", config_file)
    
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

            # create the output directory if nonexistant
            try:
                os.makedirs(parameters['output_dir'])
            except OSError:
                # directory already exists..
                if len(os.listdir(parameters['output_dir'])) != 0:
                    ans = None
                    while ans not in ["y", "n"]:
                        ans = raw_input("Writing to a non-empty directory. Files may be overwritten. Are you sure? [y/n] : ")
                    if ans == "n":
                        exit(1)

            # add a log file in the output dir
            log_fname = os.path.join(parameters['output_dir'], "pijp.log")
            hdlr = logging.FileHandler(log_fname)
            hdlr.setFormatter(log_formatter)
            logger.addHandler(hdlr)
            
            logger.info("========== writing log to file at %s ==========", log_fname)
            logger.info("Section : %s",section)
            logger.info("parameters : %s", json.dumps(parameters))
            # run the command
            segment(**parameters)

            logger.info("===== closing log")
            # remove the log file handler:
            logger.removeHandler(hdlr)


    logger.info("===== Successfuly finished pijpleiding")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= __doc__, formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()
    main(args.config_file)
