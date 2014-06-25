#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os
print("track name={}".format(sys.argv[1]))
for line in open(sys.argv[1]):
    if not line.startswith("chrom"):
        line = line.split()
        print("{}\t{}\t{}\t{}".format(line[0],line[1],line[2],line[5]))