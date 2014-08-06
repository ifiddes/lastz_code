import sys, os

for line in sys.stdin:
	line = line.split()
	names = set(line[3].split(","))
	line[3] = ",".join(sorted(names))
	print "\t".join(line)
