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
    parser.add_argument("--chain", action="store_true", default=True, help="Should lastz chain? (default: True)")
    parser.add_argument("--fasta", required=True, type=argparse.FileType("r"), help="fasta file of genes to be aligned. Will be split into new single-entry files.")
    parser.add_argument("--fastdir", default="ZNF", help="Directory fasta will be split into. default is 'ZNF'")
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


def split_fasta(file_handle, fastdir):
    filenames = list()
    for name, com, seq in read_fasta(file_handle):
        filename = fastdir + "/" + name + ".fasta"
        filenames.append(filename)
        outf = open(filename, "w")
        outf.write(">{}\n{}\n".format(name, seq))
        outf.close()
    return filenames


def call_lastz_32(call):
    """
    Generic function that runs lastz
    And a path pointing to a 1-entry fasta file containing a query chr_seq
    """
    query_path, chr_seq = call
    cmd = ["lastz_32", "--format=general:name1,zstart1,end1,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx-", "--ambiguous=n", "/dev/stdin/", query_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate(chr_seq)
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return cmdout


def call_lastz_32_chain(call):
    """
    Generic function that runs lastz with --chain on a string
    And a path pointing to a 1-entry fasta file containing a query chr_seq
    """
    query_path, chr_seq = call
    cmd = ["lastz_32", "--format=general:name1,zstart1,end1,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx-", "--chain", "--ambiguous=n", "/dev/stdin", query_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate(chr_seq)
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return cmdout


def call_generator(queries, chr_seq):
    """
    Generator to iterate over when doing multiprocessing pool
    """
    for query_path in queries:
        yield query_path, chr_seq


def main(args):
    args = parse_args(args)
    if args.fastdir not in os.listdir('.'):
        p = subprocess.Popen(["mkdir", args.fastdir]).communicate()
    else:
        print("Warning: {} already exists. Continuing with its contents.".format(args.fastdir))
    filenames = split_fasta(args.fasta, args.fastdir)
    result_dict = dict()
    pool = Pool(processes = args.cores)
    for chrom, com, seq in read_fasta(args.target):
        call_gen = call_generator(filenames, seq)
        if args.chain is True:
            result_dict[chrom] = pool.map(call_lastz_32_chain, call_gen)
        else:
            result_dict[chrom] = pool.map(call_lastz_32, call_gen)
    pool.close()
    pool.join()
    for chrom in result_dict:
        print("#{}".format(chrom))
        for result in result_dict[chrom]:
            print(result)



if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
