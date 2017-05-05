"""
Simple script to split a fasta into all individual sequences, files named based on fasta name
"""
import sys
from pyfasta import Fasta

if __name__ == '__main__':
	for name, seq in Fasta(sys.argv[1]).iteritems():
		assert ' ' not in name
		with open(name + '.fa', 'w') as outf:
			outf.write('>{}\n{}\n'.format(name, seq))
