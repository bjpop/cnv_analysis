'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : MIT
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Merge CNVs from different samples together using the collected components
of their overlap graph.
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import networkx as nx
import csv
from collections import namedtuple
from itertools import combinations


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "merge_cnvs"
DEFAULT_OVERLAP = 0.7


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Read one or more FASTA files, compute simple stats for each file'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '--overlap',
        metavar='N',
        type=float,
        default=DEFAULT_OVERLAP,
        help='CNV overlap requirement for equality (default {})'.format(
            DEFAULT_OVERLAP))
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('cnv_file',
                        metavar='CNV_FILE',
                        type=str,
                        help='Input CNV file')
    return parser.parse_args()


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


CNV = namedtuple('CNV', ['chrom', 'start', 'end', 'copynumber', 'genes'])


def collect_cnvs(cnv_filename):
    families = {}
    with open(cnv_filename) as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            this_family = row['master_sample_sheet_FAMILY_ID']
            if this_family not in families:
                families[this_family] = {}
            chroms = families[this_family]
            this_chrom = row['chr']
            this_start = int(row['coord_start'])
            this_end = int(row['coord_end'])
            this_copynumber = int(row['copy_number'])
            this_genes = tuple(row['gene_symbols_in_cnv'].split(';'))
            this_cnv = CNV(this_chrom, this_start, this_end, this_copynumber, this_genes)
            if this_chrom not in chroms:
                chroms[this_chrom] = set()
            chroms[this_chrom].add(this_cnv)
    return families 


def overlap(cnv1, cnv2, min_overlap):
    overlap_start = max(cnv1.start, cnv2.start)
    overlap_end = min(cnv1.end, cnv2.end)
    if overlap_start < overlap_end:
        overlap_size = float((overlap_end - overlap_start) + 1)
        cnv1_size = (cnv1.end - cnv1.start) + 1
        cnv2_size = (cnv2.end - cnv2.start) + 1
        cnv1_overlap = overlap_size / cnv1_size
        cnv2_overlap = overlap_size / cnv2_size
        return cnv1_overlap >= min_overlap or cnv2_overlap >= min_overlap
    return False


def merge_genes(gene_sets):
    result = set()
    for this_group in gene_sets:
        result = result.union(set(this_group)) 
    return tuple(result) 


def merge_cnvs(cnvs):
    assert(len(cnvs) > 0)
    start = min([c.start for c in cnvs])
    end = max([c.end for c in cnvs])
    genes = merge_genes([c.genes for c in cnvs]) 
    for example_cnv in cnvs:
        break
    chrom = example_cnv.chrom
    copynumber = example_cnv.copynumber
    return CNV(chrom, start, end, copynumber, genes)


def group_cnvs(family_cnvs, min_overlap):
    for family in family_cnvs:
        chroms = family_cnvs[family]
        for chrom in chroms:
            chrom_cnvs = chroms[chrom]
            overlap_graph = nx.Graph()
            overlap_graph.add_nodes_from(chrom_cnvs)
            for cnv1, cnv2 in combinations(chrom_cnvs, 2):
                if cnv1.copynumber == cnv2.copynumber and overlap(cnv1, cnv2, min_overlap):
                    overlap_graph.add_edge(cnv1, cnv2)
            for component in nx.connected_components(overlap_graph):
                if len(component) > 0:
                    merged = merge_cnvs(component)
                    yield((family, merged, component))


OUTPUT_HEADER = "\t".join(['chr', 'start', 'end', 'family', 'copy_number', 'genes'])

                
def group_and_print(family_cnvs, min_overlap):
    print(OUTPUT_HEADER)
    for family_id, cnv, component in group_cnvs(family_cnvs, min_overlap):
        gene_string = ';'.join(cnv.genes)
        print("\t".join([cnv.chrom, str(cnv.start), str(cnv.end), family_id, str(cnv.copynumber), gene_string]))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    family_cnvs = collect_cnvs(options.cnv_file)
    group_and_print(family_cnvs, options.overlap)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
