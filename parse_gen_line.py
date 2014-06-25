from __future__ import division, print_function
import sys, os
header=["name1","ztart1","end1","name2","zstart2+","end2+","size2","identity","coverage","length1","length2","nmatch","nmismatch","cigarx-","size_ratio"]
print(header)
for line in open(sys.argv[1]):
	if not line.startswith("#"):
		line = line.split()
		size1 = int(line[11].split('/')[0])
		size2 = int(line[12].split('/')[0])
		size_ratio = size1 / size2
		if float(line[8].rstrip('%')) >= 99.0 and float(line[10].rstrip('%')) >= 95.0 and ( size_ratio >= 0.9 or size_ratio <= 1.1):
			line.append(str(size_ratio))
			print("\t".join(line))

