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

logger = getLogger('pijp.bc_demultiplex')
debug, info = logger.debug, logger.info

def main(bc_index_file, sample_sheet, input_files, stats_file, output_dir, min_bc_quality, umi_length=0, bc_length=8, cut_length=35):
    """ this is the main function of this module. Does the splitting 
        and calls any other function.
    """
    cut_length = int(cut_length)
    bc_dict = create_bc_dict(bc_index_file)
    sample_dict = create_sample_dict(sample_sheet)
    files_dict = create_output_files(sample_dict, output_dir)

    try:
        sample_counter = Counter()
        for fastq_file in input_files:
            raw_fastq_dict = sanity_check_raw(fastq_file)
            il_barcode = raw_fastq_dict["index"]
            lane = raw_fastq_dict["lane"]
            flocell = raw_fastq_dict["flocell"]
            # run the splitter on this file, and collect the counts
            sample_counter += bc_split(bc_dict, sample_dict, files_dict, min_bc_quality, lane, il_barcode, flocell, fastq_file, int(umi_length), int(bc_length), cut_length)
    
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

def sanity_check_raw(fastq_file):
    """ derive lane and il_barcode from filename """
    raw_fastq_dict = dict()
    split_name = os.path.basename(fastq_file).split("_")
    raw_fastq_dict["index"] = split_name[0]
    raw_fastq_dict["barcode"] = split_name[1]
    raw_fastq_dict["lane"] = split_name[2]
    raw_fastq_dict["strand"] = split_name[3]
    with open(fastq_file, 'r') as f:
       first_line = f.readline()
    raw_fastq_dict["flocell"] = first_line.split(":")[2]
    assert ("index" in raw_fastq_dict["index"].lower()), "File name may not contain an index. Aborting (%r)" % fastq_file
    assert ("l" in raw_fastq_dict["lane"].lower()), "File name may not contain a lane. Aborting (%r)" % fastq_file
    assert ("_R1" in fastq_file), "File name does not contain R1. Aborting (%r)" % fastq_file
    r2_file = fastq_file.replace("_R1", "_R2")
    assert (os.path.isfile(r2_file)), "Could not find R2 (%r)" % r2_file
    return raw_fastq_dict

def get_sample(sample_dict, bc_dict, read, lane, il_barcode, flocell, umi_strt, umi_end, bc_strt, bc_end):
    barcode = str(read.seq[bc_strt:bc_end])
    cel_bc_id = bc_dict.get(barcode, None)
    key = (flocell, lane, il_barcode, cel_bc_id)
    return sample_dict.get(key ,None)

def filter_samples(sample_dict, lane, il_barcode, flocell):
    for sample in sample_dict:
        if not (sample.flocell==flocell and sample.lane==lane and sample.il_barcode==il_barcode):
            sample_dict.pop(sample)
    return sample_dict

def bc_split(bc_dict, sample_dict, files_dict, min_bc_quality, lane, il_barcode, flocell, r1_file, umi_length, bc_length, cut_length):
    """ Splits a fastq files according to barcode """
    sample_counter = Counter()
    umibc = umi_length + bc_length
    freader1 = FastqReader(r1_file)
    r1 = freader1
    r2_file = r1_file.replace("_R1", "_R2")
    r2 = FastqReader(r2_file)
    sample_dict = filter_samples(sample_dict, lane, il_barcode, flocell)
    if (len(sample_dict)>0):
      logger.info("splitting file %s", r1_file)
      for n, (read1, read2) in enumerate(izip(r1,r2)):
        
        # validate reads are the same
        assert (read1.name.split()[0] == read2.name.split()[0]), "Reads have different ids. Aborting."
        
        # check minimal length
        if len(read1.qual) < umibc:
            sample_counter['unqualified'] +=1
            ### skip to next iteration of loop!
            continue
        # check quality:
        quals = read1.qual[:(umibc)]
        if min(quals) >= int(min_bc_quality):
            ### trim read to cut_length
            if len(read2)>cut_length:
                read2 = read2[0:cut_length]
            # find and split
            umi_strt = 0
            umi_end = umi_length
            bc_strt = umi_length
            bc_end = bc_length+umi_length
            sample = get_sample(sample_dict, bc_dict, read1, lane, il_barcode, flocell, umi_strt, umi_end, bc_strt, bc_end)
            if (sample is not None):
                fh = files_dict[sample]
                ### ADD UMIs to the read name
                if umi_length == 0 :
                    name = read2.name
                else: 
                    ###  According to the SAM format specs, spaces are not allowed.
                    ###  Bowtie only keeps the first part, so the umi must be there.
                    name = read2.name.split()[0] + ':UMI:%s:' % read1.seq[umi_strt:umi_end]
                read2.name = name
                read = read2
                read.write_to_fastq_file(fh)
                
                sample_counter[sample] += 1
            else:
                bc = read1.seq[bc_strt:bc_end]
                fh1 = files_dict['unknown_bc_R1']
                read1.write_to_fastq_file( fh1)
                fh2 = files_dict['unknown_bc_R2']
                read2.write_to_fastq_file( fh2)

                sample_counter['undetermined'] += 1
        else:
            sample_counter['unqualified'] +=1
    else:
        logger.info("skipping file %s", r1_file)
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
    parser.add_argument('--umi-length', metavar='N', type=int, default=0,
                        help='Length of UMI (default=0, e.g. no UMI)')
    parser.add_argument('--bc-length', metavar='N', type=int, default=8,
                        help='Length of CELSeq barcode (default=8)')
    parser.add_argument('--cut-length', metavar='N', type=int, default=35,
                        help='Length of mapped read (default=35)')
    args = parser.parse_args()
    main(args.bc_index, args.sample_sheet, args.fastq_files, stats_file=args.stats_file,
         output_dir=args.out_dir, min_bc_quality=args.min_bc_quality,
         umi_length=args.umi_length, bc_length=args.bc_length, cut_length=args.cut_length)

