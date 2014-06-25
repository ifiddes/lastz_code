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

def parse_line(line):
    """
    order: znf_start, znf_end, chrom, start, end, score, identity, size_ratio, alignment_length
    """
    return [line[2], line[3], line[5], line[6], line[7], line[9]]


def main(args):
    args = parse_args(args)
    alignments = list()
    for line in args.input:
        if line.startswith("#"):
            gene = line.rstrip().lstrip("#")
            this_gene = [gene, []]
        else:
            this_gene[1].append(parse_line(line))
        alignments.append(this_gene)



if __name__ == "__main__":
    sys.exit(main(sys.argv))