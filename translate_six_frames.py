from Bio import SeqIO
import sys, os
record = SeqIO.parse(sys.argv[1], "fasta").next()
results = {}
for i in range(3):
    results["-"+str(i)] = f.reverse_complement().seq[i:].translate()
    results["+"+str(i)] = f.seq[i:].translate()

with open(sys.argv[2], "w") as outf:
    for name, seq in results.iteritems():
        outf.write(">{}_{}\n".format(record.name, name))
