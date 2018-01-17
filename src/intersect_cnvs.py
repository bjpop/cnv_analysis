'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : MIT
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Intersect CNVS
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv
from collections import namedtuple
from bx.intervals.intersection import Interval, IntervalTree
import os
from pathlib import Path 


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "intersect_cnvs"


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
    description = 'Simple case-control analysis for CNVs in families'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--genes',
                        metavar='GENES_FILE',
                        type=str,
                        help='genes coordinates')
    parser.add_argument('--cnvs',
                        metavar='CNV FILE',
                        type=str,
                        help='CNVs')
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


def read_genes(filename):
    chroms = {}
    with open(filename) as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            try:
                chrom = row['chromosome']
                start = int(row['GRCh37 start'])
                end = int(row['GRCh37 end'])
                symbol = row['symbol']
                tier = int(row['tier'])
                if chrom not in chroms:
                    chroms[chrom] = IntervalTree()
                chroms[chrom].insert(start, end, (symbol, tier))
            except:
                pass
    return chroms


def filter_cnvs(genes, filename):
    with open(filename) as file:
        reader = csv.DictReader(file, delimiter="\t")
        fieldnames = reader.fieldnames + ["tier 1 genes", "tier 2 genes", "tier 3 genes"]
        writer = csv.DictWriter(sys.stdout, fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            chrom = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            if chrom in genes:
                matches = genes[chrom].find(start, end)
                if len(matches) > 0:
                    tier1s = [ s for (s, t) in matches if t == 1]
                    tier2s = [ s for (s, t) in matches if t == 2]
                    tier3s = [ s for (s, t) in matches if t == 3]
                    row["tier 1 genes"] = ";".join(tier1s)
                    row["tier 2 genes"] = ";".join(tier2s)
                    row["tier 3 genes"] = ";".join(tier3s)
                    writer.writerow(row)



def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    genes = read_genes(options.genes)
    filter_cnvs(genes, options.cnvs)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
