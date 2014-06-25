#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, pipes, itertools, argparse
from multiprocessing import Pool

def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser(prog="pipeline",description=__doc__, usage= \
            "pipeline [args]")
    parser.add_argument("--target", required=True, type=argparse.FileType("r"), help="path to target assembly (fasta)")
    parser.add_argument("--cores", default=1, type=int, help="Number of cores (default: 1)")
    parser.add_argument("--fastadir", default="ZNF", help="Directory containing single-entry fasta files. (default: ZNF)")
    args = parser.parse_args()
    return args


def read_fasta(file_handle):
    """Generator to yield tuple representing one entry in a fasta file
    tuple format: (<id>,<comments>,<sequence>)
    Source: BioPython source code"""
    name, comments, seq = None, None, []
    for line in file_handle:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, comments, ''.join(seq))
            line = line[1:].split()
            name, comments, seq = line[0], line[1:], []
        else:
            line = ''.join([x for x in line if not x.isdigit() and not x.isspace()])
            seq.append(line)
    if name:
        yield (name, comments, ''.join(seq))


def call_lastz_32(call):
    """
    Generic function that runs lastz with --chain on a fasta entry from read_fasta
    And a path pointing to a 1-entry fasta file containing a query
    """
    query_path, (name, com, seq) = call
    cmd = ["lastz_32", "--format=general:name1,zstart1,end1,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx-", "--chain", "--ambiguous=n", query_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate(seq)
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return (name, cmdout)


def main(args):
    args = parse_args(args)
    query_files = os.listdir(args.fastadir)
    target_seqs = [x for x in read_fasta(args.target)]
    call_iterable = itertools.product(query_files, target_seqs)
    pool = Pool(processes = args.cores)
    results = pool.map(call_lastz_32, call_iterable)
    pool.close()
    pool.join()
    for name, align in results:
        sys.stdout.write("#{}\n".format(name))
        sys.stdout.write(align)




if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
