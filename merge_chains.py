#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os


def build_gene_dict(filename):
    gene_dict = dict()
    for line in open(filename):
        if line.startswith("#") and "name" not in line:
            gene = line.lstrip("#").rstrip()
            gene_dict[gene] = list()
        elif "name" not in line:
            line = line.split()
            chrom = line[0]
            chrom_start = int(line[1])
            chrom_end = int(line[2])
            identity = float(line[8].rstrip("%"))
            coverage = float(line[10].rstrip("%"))
            gene_dict[gene].append([chrom,chrom_start,chrom_end,identity,coverage])
    return gene_dict


def merge_chains(sorted_chains):
    merged = list()
    for i in xrange(1,len(sorted_chains)):
        prev_chain, next_chain = sorted_chains[i-1], sorted_chains[i]
        if prev_chain[0] == next_chain[0] and prev_chain[2] >= next_chain[1]:
            start, end = prev_chain[1], next_chain[2]
            chrom = prev_chain[0]
            prev_size, next_size = prev_chain[2] - prev_chain[1], next_chain[2] - next_chain[1]
            identity = ( prev_size * prev_chain[3] + next_size * next_chain[3] ) / ( prev_size + next_size )
            coverage = (prev_size * prev_chain[4] + next_size * next_chain[4] ) / (prev_size + next_size )
            identity, coverage = round(identity, 3), round(coverage, 3)
            merged.append([chrom, start, end, identity, coverage])
        elif i - 1 == len(sorted_chains): #if we don't merge the last chain, need to add both
            merged.append(prev_chain)
            merged.append(next_chain)
        else:
            merged.append(prev_chain)
    return merged

print("chrom\tstart\tend\tidentity\tcoverage\tgene")
gene_dict = build_gene_dict(sys.argv[1])
for gene in gene_dict:
    sorted_chains = sorted(gene_dict[gene], key = lambda x: (x[0], x[1]))
    merged = merge_chains(sorted_chains)
    highest_cov_chain = sorted(merged, key = lambda x: x[4])[-1]
    print("\t".join(map(str, highest_cov_chain + [gene])))
    #longest_chain = sorted(merged, key = lambda x: int(x[2]) - int(x[1]))[-1]
    #print("\t".join(map(str, longest_chain + [gene])))
