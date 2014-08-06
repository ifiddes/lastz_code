#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, pipes, itertools, argparse
from multiprocessing import Pool

def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser(prog="pipeline",description=__doc__, usage= \
            "pipeline [args]")
    parser.add_argument("--cores", required=True, type=int, help="number of cores")
    parser.add_argument("--file", required=True, type=argparse.FileType("r"), help="text file - list of file/file paths to be passed to shell script")
    parser.add_argument("--script", required=True, type=str, help="shell script being wrapped")
    args = parser.parse_args()
    return args

def call_chain(f):
    f, script = f
    cmd = ["sh", script, f]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate()
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)

def main(args):
    args = parse_args(args)
    pool = Pool(processes = args.cores)
    f = [x.strip() for x in args.file]
    f = [[x, args.script] for x in f]
    pool.map(call_chain, f)
    pool.close()
    pool.join()


if __name__ == "__main__":
    sys.exit(main(sys.argv))