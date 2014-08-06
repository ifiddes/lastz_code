#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os

config = """[clustalw2]
path = /cluster/home/ifiddes/sate/bin/clustalw2

[mafft]
path = /cluster/home/ifiddes/sate/bin/mafft

[muscle]
path = /cluster/home/ifiddes/sate/bin/muscle

[opal]
path = /cluster/home/ifiddes/sate/bin/opal.jar

[prank]
path = /cluster/home/ifiddes/sate/bin/prank

[raxml]
path = /cluster/home/ifiddes/sate/bin/raxml

[fasttree]
path = /cluster/home/ifiddes/sate/bin/fasttree

[commandline]
datatype = dna
input = {0}
treefile = {1}
job = satejob

[sate]
time_limit = 86400
max_subproblem_size = 0.2
aligner = mafft
merger = muscle
tree_estimator = raxml
output_directory = {2}
num_cpus = {3}
break_strategy = centroid
"""

fastas = [x.rstrip() for x in open(sys.argv[1])]

os.mkdir("sate_config_files")

for f in fastas:
    name = f.split("/")[-1].rstrip(".fasta")
    c = config.format(f, name + ".tre", "trees", 20)
    outf = open("sate_config_files/" + name + ".config", "w")
    outf.write(c)
    outf.close()