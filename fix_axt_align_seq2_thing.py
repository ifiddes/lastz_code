#!/usr/bin/env python2.7

import sys, os

f = open(sys.argv[1])


for line in f:
    if line.startswith("#") and not line.startswith("# ") and len(line) > 5:
        gene = line.lstrip("#").rstrip()
    elif line.startswith("#"):
        continue
    elif 'seq2' in line:
        line = line.replace("seq2", gene)
        sys.stdout.write(line)
    else:
        sys.stdout.write(line)

