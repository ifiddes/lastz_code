#!/usr/bin/env python

"""
Hacky file to make a .sizes file from a fasta
"""
import sys, os

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


fa = open(sys.argv[1])
for name, com, seq in read_fasta(fa):
    print "\t".join([name, str(len(seq))])

