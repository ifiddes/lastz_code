#!/usr/bin/env python2.7

import sys, os
from itertools import izip

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

records = [x for x in read_fasta(open(sys.argv[1]))]

beds = [x.split() for x in open(sys.argv[2])][1:]

for record, bed in izip(records, beds):
    name, com, seq = record
    name = ">{} {}".format(name.split("(")[0], bed[3])
    print(name)
    print(seq)
