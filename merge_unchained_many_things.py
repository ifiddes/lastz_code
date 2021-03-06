#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, argparse, itertools

def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--overlaps", type=int, required=True, nargs="+")
    parser.add_argument("--minsizes", type=int, required=True, nargs="+")
    parser.add_argument("--files", type=argparse.FileType("r"), nargs="+")
    args = parser.parse_args()
    return args


def build_gene_dict(handle):
    gene_dict = dict()
    for line in handle:
        if line.startswith("#") and "name" not in line:
            gene = line.lstrip("#").rstrip()
            gene_dict[gene] = list()
        elif "name" not in line:
            line = line.split("\t")
            chrom = line[0]
            chrom_start = int(line[1])
            chrom_end = int(line[2])
            znf_start = int(line[4])
            znf_end = int(line[5])
            identity = float(line[8].rstrip("%"))
            coverage = float(line[10].rstrip("%"))
            size_ratio = int(line[11]) / int(line[12])
            cigar = line[-1]
            gene_dict[gene].append([chrom,chrom_start,chrom_end,gene,znf_start,znf_end,identity,coverage,size_ratio,cigar])
    return gene_dict

def overlap_alignments(align, overlap):
    """
    Takes a list of alignments from lastz and tries to merge based on overlap parameter
    The list of alignments should be in the format created by build_gene_dict()
    source: http://stackoverflow.com/questions/24317211/merge-overlapping-numeric-ranges-into-continuous-ranges/24317553#24317553
    """
    # create a list of starts and ends
    stends = [ (a[0], 1) for a in align ]
    stends += [ (a[1] + overlap, -1) for a in align ]
    stends.sort(key=lambda x: x[0])
    # now we should have a list of starts and ends ordered by position,
    # e.g. if the ranges are 5..10, 8..15, and 12..13, we have
    # (5,1), (8,1), (10,-1), (12,1), (13,-1), (15,-1)
    # next, we form a cumulative sum of this
    s = 0
    cs = []
    for se in stends:
        s += se[1]
        cs.append((se[0], s))
    # this is, with the numbers above, (5,1), (8,2), (10,1), (12,2), (13,1), (15,0)
    # so, 5..8 belongs to one range, 8..10 belongs to two overlapping range,
    # 10..12 belongs to one range, etc
    # now we'll find all contiguous ranges
    # when we traverse through the list of depths (number of overlapping ranges), a new
    # range starts when the earlier number of overlapping ranges has been 0
    # a range ends when the new number of overlapping ranges is zero 
    prevdepth = 0
    start = 0
    combined = []
    for pos, depth in cs:
        if prevdepth == 0:
            start = pos
        elif depth == 0:
            combined.append((start, pos-overlap))
        prevdepth = depth
    return combined


def overlap_wrapper(gene, alignments, overlap):
    """
    Wraps around overlap_alignments(); parses gene_dict into interval list
    Returns a list that is BED format
    """
    chrom_dict = dict()
    #build a new dict mapping each chrom to a list of intervals
    for align in alignments:
        if align[0] not in chrom_dict:
            chrom_dict[align[0]] = list()
        chrom_dict[align[0]].append([align[1],align[2]])
    combined = dict()
    for chrom in chrom_dict:
        combined[chrom] = overlap_alignments(chrom_dict[chrom], overlap)
    results = list()
    for chrom in combined:
        for interval in combined[chrom]:
            results.append([chrom, interval[0], interval[1], gene])
    return results

def merge_filter_make_bed(gene_dict, overlap, minsize):
    result = list()
    for gene in gene_dict:
        merged = overlap_wrapper(gene, gene_dict[gene], overlap)
        count = 0
        for x in merged:
            if x[2] - x[1] >= minsize:
                x[-1] = x[-1] + "_" + str(count)
                result.append("\t".join(map(str, x)) + "\n")
                count += 1
    return result


def main(args):
    args = parse_args(args)
    combinations = itertools.product(args.overlaps, args.minsizes)
    for f in args.files:
        name = str(f).split("'")[1]
        gene_dict = build_gene_dict(f)
        for overlap, minsize in combinations:
            merged = merge_filter_make_bed(gene_dict, overlap, minsize)
            outf = open(name + "_overlap" + str(overlap) + "_minsize" + str(minsize) + ".bed", "w")
            for x in merged:
                outf.write(x)
            outf.close()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
