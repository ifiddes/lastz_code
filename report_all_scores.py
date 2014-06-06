#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, pipes, itertools, argparse
from multiprocessing import Pool

def build_alignment_dict(cigar_file):
    alignment_dict = dict()
    for line in cigar_file:
        if line.startswith("#"):
            gene = line.rstrip().lstrip("#")
            alignment_dict[gene] = list()
        else:
            line = line.split()
            entry = dict()
            entry['znf_start'] = line[2]
            entry['znf_end'] = line[3]
            entry['chrom'] = line[5]
            entry['start'] = line[6]
            entry['end'] = line[7]
            entry['score'] = line[9]
            entry['cigar'] = line[10:]
            entry['identity'] = calculate_identity(entry)
            entry['size_ratio'] = size_ratio(entry)
            alignment_dict[gene].append(entry)
    return alignment_dict

def calculate_identity(entry):
    """expects an entry from an entry dict (build_entry_dict)"""
    match, nonmatch = 0, 0
    cigar = entry["cigar"]
    for i in xrange(1,len(cigar),2):
        if cigar[i-1] == "M":
            match += int(cigar[i])
        else:
            nonmatch = int(cigar[i])
    return 100 * ( match / (match + nonmatch) )

def size_ratio(entry):
    """calculates ratio of aligned region in source to target. >1 = bigger in human"""
    human_size = int(entry['znf_end']) - int(entry['znf_start'])
    target_size = int(entry['end']) - int(entry['start'])
    return abs( human_size / target_size )

def parse_entry(entry):
    return [entry['score'], entry['identity'], entry['size_ratio'], entry['chrom'], entry['start'], entry['end'], entry['znf_start'], entry['znf_end'], entry['cigar']]


def main(args):
    alignments = build_alignment_dict(open(sys.argv[1]))
    header = ["gene", "score", "identity", "size_ratio", "chromosome", "start", "end", "znf_start", "znf_end", "cigar"]
    print("\t".join(header))
    for gene in alignments:
        count = 0 #give alignments unique ID
        for entry in alignments[gene]:
            result = [gene + "_" + str(count)] + parse_entry(entry)
            result[-1] = " ".join(result[-1])
            print("\t".join(map(str,result)))
            count += 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))