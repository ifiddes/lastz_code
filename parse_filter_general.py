#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, itertools, argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=argparse.FileType("r"), help="input cigar file")
    parser.add_argument("--identity", type=float, default=99, help="percent identity cutoff (default = 99)")
    parser.add_argument("--length", type=int, default=500, help="alignment length cutoff (default = 500)")
    parser.add_argument("--ratio_min", type=float, default=0.9, help="minimum size ratio (default = 0.9)")
    parser.add_argument("--ratio_max", type=float, default = 1.1, help="maximum size ratio (default = 1.1)")
    parser.add_argument("--coverage", type=float, default = 90, help="minimum coverage (default = 90)")
    parser.add_argument("--all", action="store_true", default = False, help="parse ALL chains")
    args = parser.parse_args()
    return args



def main(args):
    print("gene\tchromsome\tchrom_start\tchrom_end\tgene_start\tgene_end\tznf_size\tidentity\tcoverage\tlen_chrom\tlen_gene\tsize_ratio\tmatch\tmismatch\tcigar")
    args = parse_args(args)
    for line in args.input:
    	if line.startswith("#") and "name" not in line:
    			gene = line.lstrip("#").rstrip()
    	elif "name" not in line:
    		line = line.split()
    		chrom = line[0]
    		chrom_start = line[1]
    		chrom_end = line[2]
    		gene_start = line[4]
		gene_end = line[5]
		size2 = line[7]
		identity = float(line[8].rstrip("%"))
		coverage = float(line[10].rstrip("%"))
		len1 = int(line[11])
		len2 = int(line[12])
		size_ratio = len1 / len2
		match = line[13]
		mismatch = line[14]
		cigar = line[15]
		if identity >= args.identity and size_ratio >= args.ratio_min and size_ratio <= args.ratio_max and len1 >= args.length and len2 >= args.length and coverage >= args.coverage or args.all:
			print("\t".join(map(str,[gene,chrom,chrom_start,chrom_end,gene_start,gene_end,size2,identity,coverage,len1,len2,size_ratio,match,mismatch,cigar])))             


if __name__ == "__main__":
    sys.exit(main(sys.argv))
