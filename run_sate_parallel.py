#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess
from multiprocessing import Pool


def parallel(f):
    name = f.split("/")[-1].rstrip(".fasta")
    treepath = "trees/" + name + "/"
    os.mkdir(treepath)
    call = ["python", "/cluster/home/ifiddes/sate/run_sate.py", "-i", f, "--auto", "-j", name, "--exportconfig=" + treepath + ".config", "-o", treepath]
    subprocess.Popen(call).communicate()
    call = ["python", "/cluster/home/ifiddes/sate/run_sate.py", treepath + ".config"]
    subprocess.Popen(call).communicate()

def main(args):
    fastas = [x.rstrip() for x in open(args[1])]
    os.mkdir("trees")
    pool = Pool(processes=20)
    pool.map(parallel,fastas)
    pool.close()
    pool.join()

if __name__ == "__main__":
    sys.exit(main(sys.argv))