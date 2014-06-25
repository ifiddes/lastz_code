#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, pipes, itertools, argparse
from multiprocessing import Pool

def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser(prog="pipeline",description=__doc__, usage= \
            "pipeline [args]")
    parser.add_argument("--target", required=True, type=str, help="path to target assembly (2bit)")
    parser.add_argument("--cores", required=True, type=int, help="Number of cores")
    parser.add_argument("--output", type=argparse.FileType("w"), help="File to write CIGAR strings to")
    parser.add_argument("--chain", action="store_true", help="Should lastz chain?")
    parser.add_argument("--fasta", required=True, type=argparse.FileType("r"), help="fasta file of genes to be aligned. Mutually exclusive with --names.")
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


def build_seq_dict_fasta(fasta_handle):
    gene_sequences = dict()
    for name, com, seq in read_fasta(fasta_handle):
        gene_sequences[name] = seq
    return gene_sequences


def call_lastz_32(call):
    """
    Generic function that runs lastz on a string
    And a path pointing to a 1-entry fasta file (usually a whole genome)
    Returns the input in the lastz mapping- format
    """
    target_path, sequence, gene, chain = call
    if chain is False:
        cmd = ["lastz_32", "--format=general:name1,zstart1,end1,name2,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx-", "--ambiguous=N", target_path + "[multiple]"]
    else:
        cmd = ["lastz_32", "--format=general:name1,zstart1,end1,name2,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx-", "--chain", "--ambiguous=N", target_path + "[multiple]"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate(sequence)
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return (gene, cmdout)


def parallel_lastz(gene_sequences, target, cores, chain):
    call_iterable = list()
    for gene in gene_sequences:
        call_iterable.append((target, gene_sequences[gene], gene, chain))
    pool = Pool(processes = cores)
    alignments = pool.map(call_lastz_32, call_iterable)
    pool.close()
    pool.join()
    return alignments


def main(args):
    args = parse_args(args)
    gene_sequences = build_seq_dict_fasta(args.fasta)
    alignments = parallel_lastz(gene_sequences, args.target, args.cores, args.chain)
    if args.output is not None:
        for gene, alignment in alignments:
            args.output.write('#'+gene+'\n')
            args.output.write(alignment)
    else:
        for gene,alignment in alignments:
            sys.stdout.write('#'+gene+'\n')
            sys.stdout.write(alignment)



if __name__ == "__main__":
    sys.exit(main(sys.argv))
