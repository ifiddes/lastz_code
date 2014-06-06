#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, pipes, itertools, argparse
from multiprocessing import Pool


def read_fasta(file_handle):
    """Generator to yield tuple representing one entry in a fasta file
    tuple format: (<id>,<comments>,<sequence>)
    Source: BioPython source code"""
    name, comments, seq = None, None, []
    for line in file_handle:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, comments, ''.join(seq))
            line = line[1:].split()
            name, comments, seq = line[0], line[1:], []
        else:
            line = ''.join([x for x in line if not x.isdigit() and not x.isspace()])
            seq.append(line)
    if name:
        yield (name, comments, ''.join(seq))


def build_fasta_sizes(fasta_handle):
    size_dict = dict()
    for name, com, seq in read_fasta(fasta_handle):
        size_dict[name] = len(seq)
    return size_dict


def build_alignment_dict(cigar_file, size_dict):
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
            entry['alignment_length'] = calculate_length(entry, size_dict[gene])
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


def calculate_length(entry, size):
    """calculates coverage from an entry"""
    return str(abs(int(entry['znf_start']) - int(entry['znf_end'])))


def size_ratio(entry):
    """calculates ratio of aligned region in source to target. >1 = bigger in human"""
    human_size = int(entry['znf_end']) - int(entry['znf_start'])
    target_size = int(entry['end']) - int(entry['start'])
    return abs( human_size / target_size )


def pick_high_score(entry_list):
    sort = sorted(entry_list, key = lambda x: x['score'])
    best = sort[-1]
    return [best['score'], best['identity'], best['alignment_length'], len(sort), best['size_ratio'], best['chrom'], best['start'], best['end'], best['znf_start'], best['znf_end'], best['cigar']]


def main(args):
    fasta_sizes = build_fasta_sizes(open(sys.argv[2]))
    alignments = build_alignment_dict(open(sys.argv[1]), fasta_sizes)
    header = ["gene", "score", "identity", "alignment_length", "num_align", "size_ratio", "chromosome", "start", "end", "znf_start", "znf_end", "cigar"]
    print("\t".join(header))
    for gene in alignments:
        result = [gene] + pick_high_score(alignments[gene])
        result[-1] = " ".join(result[-1])
        print("\t".join(map(str,result)))


if __name__ == "__main__":
    sys.exit(main(sys.argv))