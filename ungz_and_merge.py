""" gunzip and merge input files to one output file.  uses plenty of memory."""

import gzip, multiprocessing
import os


def gunzip_file(input_file):
    try:
        with gzip.open(input_file, 'rb') as gzh:
            out = gzh.read()
        return out
    except KeyboardInterrupt:
        # keyboard interrupt is best managed at parent process.
        return



def main(input_files, output_dir, output_name, processes = 1):
    pool = multiprocessing.Pool(processes =processes)
    result = pool.map(gunzip_file, input_files)
    with open(os.path.join(output_dir, output_name),"wb") as fh :
        fh.write(results)
