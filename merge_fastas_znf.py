#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os
"""
Merges fasta files after the running of the net chain pipeline (net_chain_no_recip.sh)
Arguments: a text file listing all the genomes
Expects directory structure to be $GENOME/$GENOME.fasta

Outputs a fasta for each individual ZNF (maps genome/znf -> znf/genome)

Also outputs a combined fasta with every sequence pulled down
"""

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


genomes = [x.rstrip() for x in open(sys.argv[1])]

fasta_dict = dict()
combined = list()

for genome in genomes:
	gdir = os.path.join(genome, genome + ".fasta")
	for name, znf, seq in read_fasta(open(gdir)):
		if znf[0] not in fasta_dict:
			fasta_dict[znf[0]] = list()
		fasta_dict[znf[0]].append((name, genome, seq))
        

for znf in fasta_dict.iterkeys():
    for entry in fasta_dict[znf]:
        name, genome, seq = entry
        combined.append(">{}_{}_{}\n{}\n".format(name, genome, znf, seq))


os.mkdir("znf_fastas_for_trees")
for znf in fasta_dict.iterkeys():
	outf = open(os.path.join("znf_fastas_for_trees", znf + ".fasta"), "w")
	for name, genome, seq in fasta_dict[znf]:
		outf.write(">{} {}\n{}\n".format(name, genome, seq))
	outf.close()


outf = open("everything_combined.fasta", "w")
for line in combined:
    outf.write(line)

outf.close()

