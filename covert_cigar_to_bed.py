#!/usr/bin/env python2.7

from __future__ import print_function
import sys

for line in open(sys.argv[1]):
    line = line.split()
    if not line[0].startswith("gene"):
        if float(line[2]) > 99 and int(line[3]) > 2000 and float(line[5]) > 0.9 and float(line[5]) < 1.1:
            outline = [line[6], line[7], line[8], line[0]]
            print("\t".join(outline))

