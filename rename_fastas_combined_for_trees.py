#!/usr/bin/env python
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

fastas = list()

#for f in open(sys.argv[1]):
for f in open("files"):
	f = f.rstrip()
	g = os.path.basename(f).split(".")[0]
	for name, com, seq in read_fasta(open(f)):
		name = name + "_" + g
		name = name.replace(":",";")
		fastas.append(">{}\n{}\n".format(name,seq))

outf = open("combined.renamed.fasta","w")
for x in fastas:
	outf.write(x)

