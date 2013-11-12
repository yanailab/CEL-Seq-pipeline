#!/usr/bin/python2
""" Clean up script after the previous segments. Deletes whatever files are given
in `input_files`, and empty parent directories of these files.

USERS BEWARE:  really deletes everything without confirmation!
"""

from logging import getLogger
import os

logger = getLogger("pijp.cleanup")

def main(input_files, really_delete):
    
    if not really_delete:
        logger.info("Running clean_up in 'demo' mode - not really deleting anything")

    #  Go over the supplied filenames and delete them.
    for filename in input_files:
        logger.info("Deleting %s", filename)
        if really_delete:
            try:
                os.remove(filename)
            except OSError as e:
                logger.error("Error deleting %s: %s", filename, e.strerror) 

    #  Go over the parent dirs and delete them.
    #  Recursively deletes empty parents.
    for dirname in set(os.path.dirname(x) for x in input_files):
        logger.info("Trying to delete %s and empty parent dirs", dirname)
        if really_delete:
            try:
                os.removedirs(dirname)
            except OSError:
                pass  # the dir is probabily not empty. 
