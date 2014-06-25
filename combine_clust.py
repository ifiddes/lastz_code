#!/usr/bin/env python2.7

import sys, os
resultdirs = [x for x in os.listdir('.') if x.endswith('.fasta')]
combined = list()
for rd in resultdirs:
	results = [x for x in os.listdir(rd) if x.endswith('.result')]
	for r in results:
		for line in open(rd+"/"+r):
			if not line.startswith('#'):
				combined.append(line)

header = "#name1	zstart1	end1	name2	zstart2+	end2+	size2	identity	idPct	coverage	covPct	length1	length2	nmatch	nmismatch	cigarx-"

sys.stdout.write(header)
sys.stdout.write("\n")
for line in combined:
	sys.stdout.write(line)
