#!/usr/bin/python2
""" Paired end barcode splitter

Splits one fastq file based on barcodes from another file. The script 
accepts one barcode index file, one sample sheet, and a set of fastq files.
The fastq files are assumed to be read1 (the barcode), and the script replaces
"_R1_" with "_R2_" to find the second file. If the filename ends in `gz` it
is decompressed on-the-fly (HTSeq does that).

"""

from __future__ import print_function, division

import os
import argparse
import csv
from logging import getLogger
from itertools import izip
from collections import OrderedDict, namedtuple, Counter

from HTSeq import FastqReader, SequenceWithQualities

FN_SCHEME = "{0.project}_{0.series}_sample_{0.id}.fastq"
FN_UNKNOWN = "undetermined_{0}.fastq"
CUT_LENGTH = 35


logger = getLogger('pijp.bc_demultiplex')
debug, info = logger.debug, logger.info

def main(bc_index_file, sample_sheet, input_files, stats_file, output_dir, min_bc_quality, umi_length=0, bc_length=8, kleine=False):
    """ this is the main function of this module. Does the splitting 
        and calls any other function.
    """

    bc_dict = create_bc_dict(bc_index_file)
    sample_dict = create_sample_dict(sample_sheet)
    files_dict = create_output_files(sample_dict, output_dir)

    if kleine is True:
        #small seq
        bc_split = bc_split_se
    else:
        bc_split = bc_split_pe

    try:
        sample_counter = Counter()
    
        for fastq_file in input_files:
            logger.info("splitting file %s", fastq_file)
            r1_file = fastq_file
        
            # derive lane and il_barcode from filename
            split_name = os.path.basename(fastq_file).split("_")
            il_barcode = split_name[0]
            lane = split_name[2]
    
            # run the splitter on this file, and collect the counts
            sample_counter += bc_split(bc_dict, sample_dict, files_dict, min_bc_quality, lane, il_barcode, r1_file , int(umi_length), int(bc_length))
    
        ### Create the stats file.
        total = sum(sample_counter.values())
        stats = [["# Sample_id", "reads", "precentage"]]
        
        # a sample can appear more than once in the sample dict,
        # but we do not want to use set as it is unordered. 
        samples = OrderedDict(((x, True) for x in sample_dict.values()))
        
        for sample in samples.keys():
            sample_count = sample_counter[sample]
            stats.append( [FN_SCHEME.format(sample), sample_count, 100.0*sample_count/total])
        stats.append( ["unqualified", sample_counter['unqualified'], 100.0*sample_counter['unqualified']/total])
        stats.append( ["undetermined", sample_counter['undetermined'], 100.0*sample_counter['undetermined']/total])
        stats.append( ["total", total, 100])
        with open(os.path.join(output_dir,stats_file), "w") as stats_fh:
            stats_writer = csv.writer(stats_fh, delimiter='\t')
            stats_writer.writerows(stats)
    finally:        
        for file in files_dict.values():
            file.close()

def create_output_files(sample_dict, target):

    files_dict = dict()
    for sample in set(sample_dict.values()):
        filename = os.path.join(target, FN_SCHEME.format(sample))
        files_dict[sample] = open(filename,"wb")
    filename = os.path.join(target, FN_UNKNOWN.format('R1'))
    files_dict['unknown_bc_R1'] = open(filename, "wb")
    filename = os.path.join(target, FN_UNKNOWN.format('R2'))
    files_dict['unknown_bc_R2'] = open(filename, "wb")
    return files_dict
    
def create_sample_dict(sample_sheet_file):
    """  Create a mapping from sample keys to sample infos """
    sample_dict = OrderedDict()
    # define the "data types"
    Key = namedtuple('sample_key',['flocell', 'lane', 'il_barcode', 'cel_barcode'])
    Sample_info = namedtuple('Sample', ['id', 'series', 'project'])

    with open(sample_sheet_file, 'rb') as sample_sheet_fh:
        sample_sheet_reader = csv.DictReader(sample_sheet_fh, delimiter='\t')
        for row in sample_sheet_reader:
            id = "{0:04}".format(int(row["#id"]))  #  id has an extra "#" becaues its the first field
            key = Key(row["flocell"], row["lane"], row["il_barcode"], row["cel_barcode"])
            sample_dict[key] = Sample_info(id, row["series"], row["project"])

    return sample_dict                              

def create_bc_dict(bc_index_file):
    """ create a mapping from barcode sequence to barcode id """
    bc_dict = dict()
    with open(bc_index_file, 'rb') as bc_index:
        bc_index_reader = csv.reader(bc_index, delimiter='\t')
        for row in bc_index_reader:
            if row[0].startswith("#"):
                continue
            bc_dict[row[1]] = row[0]
    return bc_dict


def get_sample(sample_dict, bc_dict, read, lane, il_barcode, umi_length, bc_length):
    barcode = str(read.seq[umi_length:(umi_length+bc_length)])
    cel_bc_id = bc_dict.get(barcode, None)
    flocell = read.name.split(":")[2]
    key = (flocell, lane, il_barcode, cel_bc_id)
    return sample_dict.get(key ,None)


def bc_split_pe(bc_dict, sample_dict, files_dict, min_bc_quality, lane, il_barcode, r1_file, umi_length, bc_length):
    """ Splits a pair of fastq files according to barcode """
    sample_counter = Counter()
    assert ("_R1_" in r1_file), "File name does not contain R1. Aborting"
    r1 = FastqReader(r1_file)
    r2_file = r1_file.replace("_R1_", "_R2_")
    r2 = FastqReader(r2_file)
    for n, (read1, read2) in enumerate(izip(r1,r2)):
        if (n % 1e5)==0:
            debug("read number {0}".format(n))
        # validate reads are the same
        assert (read1.name.split()[0] == read2.name.split()[0]), "Reads have different ids. Aborting."
        # check quality:
        quals = read1.qual[:(umi_length+bc_length)]
        if min(quals) >= int(min_bc_quality):
            # find and split
            sample = get_sample(sample_dict, bc_dict, read1, lane, il_barcode, umi_length, bc_length)
            if (sample is not None):
                fh = files_dict[sample]
                ### ADD UMIs to the read name
                read2.name += ' UMI:%s' % read1.seq[:umi_length]
                read2.write_to_fastq_file(fh)
                sample_counter[sample] += 1
            else:
                fh1 = files_dict['unknown_bc_R1']
                fh2 = files_dict['unknown_bc_R2']
                read1.write_to_fastq_file( fh1)
                read2.write_to_fastq_file( fh2)
                sample_counter['undetermined'] += 1
        else:
            sample_counter['unqualified'] +=1
    return sample_counter

def bc_split_se(bc_dict, sample_dict, files_dict, min_bc_quality, lane, il_barcode, r1_file, umi_length, bc_length):
    """ Splits a single fastq file according to barcode """
    sample_counter = Counter()
    umibc = umi_length + bc_length
    r1 = FastqReader(r1_file)
    for n, read1 in enumerate(r1):
        if (n % 1e5)==0:
            debug("read number {0}".format(n))
        # check quality:
        
        quals = read1.qual[:(umibc)]
        if min(quals) >= int(min_bc_quality):
            # find and split
            sample = get_sample(sample_dict, bc_dict, read1, lane, il_barcode, umi_length, bc_length)
            if (cel_bc_id is not None) and (sample is not None):
                fh = files_dict[sample]
                ### ADD UMIs to the read name
                name = read1.name + ' UMI:%s' % read1.seq[:umi_lengt]
                read = SequenceWithQualities(read1.seq[umibc:], name, read1.qual[umibc:])
                read1.write_to_fastq_file(fh)
                sample_counter[sample] += 1
            else:
                fh1 = files_dict['unknown_bc_R1']
                read1.write_to_fastq_file( fh1)
                sample_counter['undetermined'] += 1
        else:
            sample_counter['unqualified'] +=1
    return sample_counter

                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= __doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--min-bc-quality', metavar='N', type=int, default=10, 
                        help='Minimal quality for barcode reads (default=10)')
    parser.add_argument('--out-dir', metavar='DIRNAME', type=str, default='.',
                        help='Output directory. Defaults to current directory')
    parser.add_argument('--stats-file', metavar='STATFILE', type=str, default='stats.tab',
                        help='Statistics file name (default: stats.tab)')
    parser.add_argument('bc_index', type=str)
    parser.add_argument('sample_sheet', type=str)
    parser.add_argument('fastq_files', type=str, nargs='+')
    args = parser.parse_args()
    main(args.bc_index, args.sample_sheet, args.fastq_files, stats_file=args.stats_file,
         output_dir=args.out_dir, min_bc_quality=args.min_bc_quality)

