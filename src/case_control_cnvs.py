'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : MIT
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX


Simple case-control analysis of CNVs by family.

Compare:

    positive cases, negative cases, positive controls, negative controls

Using chi-squared test.

Does not work when the various counts for cases/contols is small (particularly zero).

Also removes duplicate variant calls for samples (same sample ID, different sentrix ID).
Some samples were tested more than once for QC. We should be careful not to count
them as the same thing. We arbitrarily choose to keep the first set of CNVs for a sample
and write the duplicates out to a file.
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
from pathlib import Path 
from scipy.stats import chi2_contingency


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "case_control_cnvs"


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
    parser.add_argument('--merged',
                        metavar='MERGED_CNV_FILE',
                        type=str,
                        help='Input merged CNV file (output of merge_cnvs')
    parser.add_argument('--all',
                        metavar='ALL_CNV_FILE',
                        type=str,
                        help='Input all CNV file containing all CNVs in all families')
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


def read_merged_cnvs(merged_cnvs_filename):
    families = {}
    with open(merged_cnvs_filename) as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            this_family = row['family']
            if this_family not in families:
                families[this_family] = {}
            chroms = families[this_family]
            this_chrom = row['chr']
            this_start = int(row['start'])
            this_end = int(row['end'])
            this_copynumber = int(row['copy_number'])
            this_genes = tuple(row['genes'].split(';'))
            this_cnv = CNV(this_chrom, this_start, this_end, this_copynumber, this_genes)
            if this_chrom not in chroms:
                chroms[this_chrom] = IntervalTree() 
            chroms[this_chrom].insert(this_start, this_end, this_cnv)
    return families 


class CasesControls(object):
    def __init__(self):
        self.cases = set()
        self.controls = set()

def read_all_cnvs(all_cnvs_filename):
    families = {}
    seen_samples = {} 
    duplicates = []
    cases_controls = {}
    with open(all_cnvs_filename) as file:
        reader = csv.DictReader(file, delimiter='\t')
        header = reader.fieldnames 
        for row in reader:
            this_family = row['master_sample_sheet_FAMILY_ID']
            this_sample_id = row['sample_id']
            if this_family not in cases_controls:
                cases_controls[this_family] = CasesControls()
            if row['ped_Affected'] == "Yes":
                cases_controls[this_family].cases.add(this_sample_id)
            else:
                cases_controls[this_family].controls.add(this_sample_id)
            this_sentrix_id = row['sentrix_id']
            if this_sample_id not in seen_samples:
                seen_samples[this_sample_id] = this_sentrix_id
            if this_sentrix_id == seen_samples[this_sample_id]:
                # not a duplicate
                if this_family not in families:
                    families[this_family] = []
                families[this_family].append(row)
            else:
                duplicates.append(row)
        return header, duplicates, families, cases_controls



def write_family_results(outdir, family_id, results):
    cnv_id = 0
    nodes = []
    unique_samples = set()
    unique_cnvs = {}
    nodes = []
    edges = []
    for sample, cnv in results:
        unique_samples.add(sample)
        if cnv not in unique_cnvs:
            unique_cnvs[cnv] = cnv_id
            cnv_id += 1
        edges.append({'data': {'source': sample.id, 'target': str(unique_cnvs[cnv])}})
    for sample in unique_samples:
        if sample.affected:
            taxon = "CASE"
        else:
            taxon = "CONTROL"
        nodes.append({'data': {'id': sample.id, 'annotation_Taxon': taxon }})
    for cnv, cnv_id in unique_cnvs.items():
        cnv_name = cnv.genes[0] 
        cnv_alias = cnv.genes
        cnv_taxon = 'CNV' + str(cnv.copynumber)
        nodes.append({'data': {'id': str(cnv_id), 'name': cnv_name, 'alias': cnv_alias, 'annotation_Taxon': cnv_taxon}})
    with open(os.path.join(outdir, family_id + '.json'), "w") as outfile:
        result = {'elements': {'nodes': nodes, 'edges': edges}}
        outfile.write(json.dumps(result)) 


SAMPLE = namedtuple('SAMPLE', ['id', 'affected'])


def intersect_cnvs(merged_cnvs, all_cnvs):
    result = {}
    for family_id, this_family_cnvs in all_cnvs.items():
        if family_id not in result:
            result[family_id] = {}
        this_family_result = result[family_id]
        this_merged_cnvs = merged_cnvs[family_id]
        for this_cnv in this_family_cnvs:
            this_affected = this_cnv['ped_Affected'] == "Yes"
            this_sample = SAMPLE(this_cnv['sample_id'], this_affected)
            this_chrom = this_cnv['chr']
            this_interval_tree = this_merged_cnvs[this_chrom]
            this_intersection = this_interval_tree.find(int(this_cnv['coord_start']), int(this_cnv['coord_end']))
            for intersecting_cnv in this_intersection:
                if intersecting_cnv not in this_family_result:
                    this_family_result[intersecting_cnv] = set()
                this_family_result[intersecting_cnv].add(this_sample)
    return result


def write_duplicates(fieldnames, duplicates, input_filename):
    input_path = Path(input_filename)
    output_filepath = input_path.with_suffix(".dups.tsv")
    with output_filepath.open("w") as output_file:
        writer = csv.DictWriter(output_file, fieldnames, delimiter="\t")
        writer.writeheader()
        for row in duplicates:
            writer.writerow(row)



def affected_str(is_affected):
    if is_affected:
        return "CASE"
    else:
        return "CONTROL"


def write_family_intersections(cnv_significances):
    header = '\t'.join(["chr", "start", "end", "family", "positive cases", "negative cases","postive controls","negative controls", "chi2", "p-value", "copynumber", "genes", "samples"])
    print(header)
    for family_id, pos_cases, neg_cases, pos_controls, neg_controls, chi2, p, chrom, start, end, copynumber, genes, samples in cnv_significances: 
        gene_string = ";".join(genes)
        sample_string = "|".join([";".join([s.id, affected_str(s.affected)]) for s in samples])
        print("\t".join([chrom, str(start), str(end), family_id, str(pos_cases), str(neg_cases), str(pos_controls), str(neg_controls), str(chi2), str(p), str(copynumber), gene_string, sample_string]))


def get_significance(cases_controls, family_intersections):
    result = []
    for family_id, cnvs in family_intersections.items():
        this_cases_controls = cases_controls[family_id]
        total_cases = len(this_cases_controls.cases)
        total_controls = len(this_cases_controls.controls)
        for this_cnv, samples in cnvs.items():
            positive_cases = len([s for s in samples if s.affected])
            positive_controls = len([s for s in samples if not s.affected])
            negative_cases = total_cases - positive_cases
            negative_controls = total_controls - positive_controls
            try:
                (chi2, p, dof, ex) = chi2_contingency([[positive_cases, positive_controls], [negative_cases, negative_controls]])
            except ValueError:
                chi2 = 0.0
                p = 1.0
            this_result = [family_id, positive_cases, negative_cases, positive_controls, negative_controls, chi2, p, this_cnv.chrom, this_cnv.start, this_cnv.end, this_cnv.copynumber, this_cnv.genes, samples]
            result.append(this_result)
    return result


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    merged_cnvs = read_merged_cnvs(options.merged)
    header, duplicates, families, cases_controls = read_all_cnvs(options.all)
    write_duplicates(header, duplicates, options.all)
    family_intersections = intersect_cnvs(merged_cnvs, families)
    cnv_significances = get_significance(cases_controls, family_intersections)
    write_family_intersections(cnv_significances)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
