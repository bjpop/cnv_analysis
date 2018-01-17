'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : MIT
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Graph the relationship between samples and CNVs.

Output is a JSON file per family, in a format suitable for Cytoscape
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import networkx as nx
import csv
from collections import namedtuple
from itertools import combinations
from bx.intervals.intersection import Interval, IntervalTree
import json
import os
import networkx as nx


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "graph_cnvs"
DEFAULT_OVERLAP = 0.7
DEFAULT_OUTPUT_DIR = "out" 


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
    description = 'Graph the relationships between samples and CNVs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('cnvs',
                        metavar='CNV_FILE',
                        type=str,
                        help='Input CNV file (output of cnv_case_control)')
    parser.add_argument('--out',
                        metavar='OUT_DIR',
                        type=str,
                        default=DEFAULT_OUTPUT_DIR,
                        help='Output directory for all family graphs in JSON format')
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


CNV = namedtuple('CNV', ['chrom', 'start', 'end', 'copynumber', 'chi2', 'p', 'penncnv_conf', 'genes'])


def read_all_cnvs(all_cnvs_filename):
    families = {}
    with open(all_cnvs_filename) as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            this_family = row['family']
            if this_family not in families:
                families[this_family] = []
            families[this_family].append(row)
    return families 


SAMPLE = namedtuple('SAMPLE', ['id', 'affected'])

def make_sample(sample_info):
    id, affected = sample_info.split(';')
    return SAMPLE(id, affected)

def write_family_results(outdir, family, cnvs):
    cnv_id = 0
    unique_samples = set()
    unique_cnvs = {}
    unique_edges = set()
    graph = nx.Graph() 
    for row in cnvs:
        this_samples = [make_sample(sample_info) for sample_info in row['samples'].split('|')]
        for s in this_samples:
            unique_samples.add(s)
        this_genes = tuple(row['genes'].split(';'))
        this_cnv = CNV(row['chr'], row['start'], row['end'], row['copynumber'], row['chi2'], row['p-value'], row['penncnv_conf'], this_genes)
        if this_cnv not in unique_cnvs:
            unique_cnvs[this_cnv] = cnv_id
            cnv_id += 1
        for sample in this_samples:
            unique_edges.add((sample.id, str(unique_cnvs[this_cnv])))
    for sample in unique_samples:
        node_id = sample.id
        if sample.affected == "CASE":
            this_name = "+"
        else:
            this_name = "-"
        graph.add_node(node_id, name=this_name, type=sample.affected)
    for cnv, cnv_id in unique_cnvs.items():
        node_id = str(cnv_id)
        node = graph.add_node(node_id, name=cnv.genes[0], genes=";".join(cnv.genes), type='CNV' + str(cnv.copynumber), p=float(cnv.p), chi2=float(cnv.chi2), penncnv_conf=float(cnv.penncnv_conf))
    for source, target in unique_edges:
        graph.add_edge(source, target)

    outfilename = os.path.join(outdir, family + ".graphml")
    nx.write_graphml(graph, outfilename, prettyprint=True, infer_numeric_types=True)


def make_output_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def make_graphs(outdir, family_cnvs):
    for family, cnvs in family_cnvs.items():
        write_family_results(outdir, family, cnvs)

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    make_output_directory(options.out)
    family_cnvs = read_all_cnvs(options.cnvs)
    make_graphs(options.out, family_cnvs)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
