#!/usr/bin/env python
"""combines individual bedgraphs into one large datatable for manual editing"""

import sys, os

files = [x for x in os.listdir(".") if x.endswith(".bedGraph") and "median" not in x and "variance" not in x]

median = "median.bedGraph"; variance = "variance.bedGraph"

from collections import OrderedDict as od

table = od()
for f in sorted(files):
    sample = os.path.basename(f).split('.')[0]
    lines = [x.split() for x in open(f)][1:]
    table[sample] = dict([[x[2], x[3]] for x in lines])

med_var = list()
med = [x.split() for x in open("median.bedGraph")][1:]
var = [x.split() for x in open("variance.bedGraph")][1:]

final = list()

for m, v in zip(med, var):
    pos = m[2]
    mval = round(float(m[3]), 4); vval = round(float(v[3]), 4)
    kvals = []
    for s in table:
        if pos in table[s]:
            kvals.append(round(float(table[s][pos]), 4))
        else:
            kvals.append("0")
    final.append([pos] + kvals + [mval, vval])

outf = open(sys.argv[1], "w")
outf.write("\t".join(["position"] + sorted(files) + ["median","variance"])); outf.write("\n")
for x in final:
    if sum(map(float, x[1:])) > 0:
        outf.write("\t".join(map(str,x))); outf.write('\n')

outf.close()
