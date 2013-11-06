#!/usr/bin/python

import ConfigParser
import argparse


import bowtie_wrapper, bc_demultiplex


def main(config_file):
    
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    sections = ("bc_demultiplex", "bowtie_wrapper", "htseq_wrapper")
    wrappers = (bc_demultiplex.main, bowtie_wrapper.main, htseq_wrapper.main)
    
    for section, wrapper in zip(sections, wrappers):

        wrapper(*dict(*conf))

