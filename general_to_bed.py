from __future__ import print_function
import sys
for line in open(sys.argv[1]):
	if line.startswith('#'):
		if 'name1' not in line:
			gene = line.lstrip('#').rstrip()
	else:
		line = line.split()
		print("\t".join([line[0],line[1],line[2],gene]))
