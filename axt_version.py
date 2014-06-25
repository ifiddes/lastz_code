#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, pipes, itertools, argparse
from multiprocessing import Pool

def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser(prog="pipeline",description=__doc__, usage= \
            "pipeline [args]")
    parser.add_argument("--target", required=True, type=str, help="target assembly. Will pull 2bit from /cluster/data/$assembly/$assembly.2bit to make separate fastas.")
    parser.add_argument("--cores", default=1, type=int, help="Number of cores (default=1")
    parser.add_argument("--output", type=argparse.FileType("w"), help="File to write CIGAR strings to. If unset, stdout is used.")
    parser.add_argument("--fasta", type=argparse.FileType("r"), help="fasta file of genes to be aligned. Mutually exclusive with --names.")
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


def call_twoBitToFa_WholeGenome(target):
    """
    Converts whole genome (target) into chromosome chunk fastas
    """
    genome_path = "/cluster/data/{}/{}.2bit".format(target, target)
    info_cmd = ["twoBitInfo", genome_path, "/dev/stdout"]
    p = subprocess.Popen(info_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmdout, cmderr = p.communicate()
    cmdout = cmdout.split()
    chromosomes = list()
    for x in cmdout:
        if 'chr' in x:
            chromosomes.append(x)
    make_dir = ["mkdir", target]
    p = subprocess.Popen(make_dir).communicate()
    for chrom in chromosomes:
        chrom_cmd = ["twoBitToFa", "-seq=" + chrom, genome_path, "{}/{}.fasta".format(target, chrom)]
        p = subprocess.Popen(chrom_cmd, stderr=subprocess.PIPE).communicate()
    return chromosomes


def call_lastz(call):
    """
    Generic function that runs lastz on a string
    And a path pointing to a 1-entry fasta file (usually a whole chromosome)
    Returns the input in the lastz mapping- format
    """
    target_path, sequence, gene = call
    cmd = ["lastz", "--format=axt+", "--ambiguous=n" , target_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate(sequence)
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return (gene, cmdout)


def parallel_lastz(gene_sequences, chromosomes, target, cores):
    call_iterable = list()
    for chrom in chromosomes:
        chrom_path = "{}/{}.fasta".format(target, chrom)
        for gene in gene_sequences:
            call_iterable.append((chrom_path, gene_sequences[gene], gene))
    pool = Pool(processes = cores)
    alignments = pool.map(call_lastz, call_iterable)
    pool.close()
    pool.join()
    return alignments


def main(args):
    args = parse_args(args)
    gene_sequences = build_seq_dict_fasta(args.fasta)
    if args.target not in os.listdir("."): #see if genome has been split before
        chromosomes = call_twoBitToFa_WholeGenome(args.target)
    else:
        #build a fasta of target genome, 1 fasta per entry
        chromosomes = [x.split(".")[0] for x in os.listdir(args.target)]
    alignments = parallel_lastz(gene_sequences, chromosomes, args.target, args.cores)
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
