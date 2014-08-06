#!/usr/bin/env python2.7
from __future__ import division, print_function
import os, sys

r = list()

l = [x.split() for x in open(sys.argv[1])]


x = l[0]
for i in xrange(1, len(l)):
    if x[3] != l[i][3]:
        r.append(x)
        x = l[i]
    elif int(l[i][2]) > int(x[2]):
        x = [x[0], x[1], l[i][2], x[3]]
    else:
        r.append(x)


genome = sys.argv[1].split(".")[0]
print('track name={0}_merged description="{0} merged best reciprocal chains"'.format(genome))
for x in r:
    print("\t".join(x))