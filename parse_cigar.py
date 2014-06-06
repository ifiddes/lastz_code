#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, itertools, argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=argparse.FileType("r"), help="input cigar file")
    parser.add_argument("--identity", type=float, help="percent identity cutoff")
    parser.add_argument("--length", type=int, help="alignment length cutoff")
    parser.add_argument("--ratio_min", type=float, help="minimum size ratio")
    parser.add_argument("--ratio_max", type=float, help="maximum size ratio")
    args = parser.parse_args()
    return args

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
            entry['alignment_length'] = calculate_length(entry)
            alignment_dict[gene].append(entry)
    return alignment_dict

def parse_entry(entry, identity, length, ratio_min, ratio_max):
    if entry['alignment_length'] >= length and entry['size_ratio'] >= ratio_min and entry['size_ratio'] <= ratio_max and entry['identity'] >= identity:
        return [entry['score'], entry['identity'], entry['alignment_length'], entry['size_ratio'], entry['chrom'], entry['start'], entry['end'], entry['znf_start'], entry['znf_end'], entry['cigar']]
    else:
        return None

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


def calculate_length(entry):
    return str(abs(int(entry['znf_start']) - int(entry['znf_end'])))


def size_ratio(entry):
    """calculates ratio of aligned region in source to target. >1 = bigger in human"""
    human_size = int(entry['znf_end']) - int(entry['znf_start'])
    target_size = int(entry['end']) - int(entry['start'])
    return abs( human_size / target_size )


def main(args):
    args = parse_args(args)
    alignments = build_alignment_dict(args.input)
    header = ["gene", "score", "identity", "alignment_length", "size_ratio", "chromosome", "start", "end", "znf_start", "znf_end", "cigar"]
    print("\t".join(header))
    for gene in alignments:
        for alignment in alignments[gene]:
            result = [gene] + parse_entry(alignment, args.identity, args.length, args.ratio_min, args.ratio_max)
            if result[-1] is not None:
                print(result)
                #result[-1] = " ".join(result)
                #print("\t".join(map(str,result)))


if __name__ == "__main__":
    sys.exit(main(sys.argv))