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

###################################################################################################
## sections are the names of sections in the config file.
## segments are the functions you run. The segments and sections MUST be ordered
## the same way.
SECTIONS = ( "bc_demultiplex", "bowtie_wrapper", "htseq_wrapper", "clean_up")
SEGMENTS = ( bc_demultiplex.main, bowtie_wrapper.main, htseq_wrapper.main, clean_up.main)
###################################################################################################

# some definitions for the loggers.

LOGFORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logger = logging.getLogger('pijp')
logging.basicConfig(level=logging.INFO, format = LOGFORMAT)
log_formatter = logging.Formatter(LOGFORMAT)

def main(config_file):
    
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    # does command do anything? called before handle is set. any point to it?
    logger.info("===== Started pijpleiding with config file : %s =====", config_file)
    
    for section, segment in zip(SECTIONS, SEGMENTS):

        if config.getboolean(section, "pipe_run"):
            parameters = dict(config.items(section))
            parameters.pop("pipe_run")  #  Remove this from the dictionary, so it will not be passed to the segment.

            ##  Read the input file list.
            ##  pipe_input_files can have multiple lines, and each line is a glob pattern,
            ##  i.e. include "*" and "?" as in regular bash.
            parameters["input_files"] = []
            for input_glob_pattern in parameters.pop("pipe_input_files").splitlines():
                gl = glob(input_glob_pattern)
                assert (len(gl) != 0 ), "input files not found for pattern "+input_glob_pattern
                parameters["input_files"] += gl
            parameters["input_files"].sort()

            ##
            create_dir(parameters['output_dir'])

            # add a log file in the output dir
            log_fname = os.path.join(parameters['output_dir'], "pijp.log")
            hdlr = logging.FileHandler(log_fname)
            hdlr.setFormatter(log_formatter)
            logger.addHandler(hdlr)
            logger.info("===== pipeline at %s =====", os.path.realpath(__file__))
            logger.info("========== writing log to file at %s ==========", log_fname)
            logger.info("Section : %s",section)
            logger.info("parameters : %s", json.dumps(parameters))


            ####  run the command ============================================
            # that's the heart of the whole pipeline
            # in example, bc_demultiplex.main(parameters_from_log_file) 
            segment(**parameters)
            #### =============================================================

            # remove the log file handler:
            logger.info("=========== closing log ===========")
            logger.removeHandler(hdlr)

    # as beofre, not seen in log file. declared after handle closed. any point to it?
    logger.info("===== Successfuly finished pijpleiding =====")




def create_dir(dirname):
    ## create the output directory if nonexistant
    try:
        os.makedirs(dirname)
    except OSError:
        # directory already exists.. check if empty or
        # perhaps has only pijp.log file. otherwise, prompt the user
        if (len(os.listdir(dirname)) != 0) and (os.listdir(dirname) != ['pijp.log']):
            ans = None
            while ans not in ["y", "n"]:
                print("Opening directory {0} ".format(dirname))
                ans = raw_input("Writing to a non-empty directory. Files may be overwritten. Are you sure? [y/n] : ")
            if ans == "n":
                logger.info("Aborted by user because of non empty dir")
                exit(1)


if __name__ == "__main__":

    ##  Parse command line options. This adds the useful --help option.
    parser = argparse.ArgumentParser(description= __doc__, formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()
    main(args.config_file)
