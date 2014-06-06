#!/usr/bin/env python2.7
from __future__ import division, print_function
import sys, os, subprocess, argparse


def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser(prog="pipeline",description=__doc__, usage= \
            "pipeline [args]")
    parser.add_argument("--assembly", required=True, type=str, help="source assembly")
    parser.add_argument("--offset", default=1000, type=int, help="Number of bases upstream and downstream to pull down (default=1000")
    parser.add_argument("--names", required=True, type=argparse.FileType("r"), help="file with gene names to be aligned. Mutually exclusive with --fasta.")
    parser.add_argument("--fasta", required=True, type=argparse.FileType("w"), help="File to write gene sequences to.")
    parser.add_argument("--failed", required=True, type=argparse.FileType("w"), help="File to write gene names that were not found to.")
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
        

def names_to_seqs(name_file, assembly, offset):
    """
    Uses hgsql to turn gene names into regions from annotation of <assembly>
    And adding <offset> bases to each side. If multiple genes are returned,
    the longest annotation is used.
    Returns a dict of regions mapping each gene to one region
    Also returns a list of genes which did not have annotation
    """
    query = "select k.chrom,k.txStart,k.txEnd,x.geneSymbol,k.exonCount,k.strand from knownGene k, kgXref x where k.name=kgID and x.geneSymbol=\"{}\"" 
    gene_regions = dict(); unknown_genes = list()
    for name in name_file:
        name = name.rstrip()
        this_gene = callHgsql(assembly, query.format(name))
        if len(this_gene) == 0:
            unknown_genes.append(name)
        else:
            this_gene = this_gene.rstrip().split("\n")
            for entry in this_gene:
                entry = entry.split()
                entry[1] = int(entry[1]) - offset
                entry[2] = int(entry[2]) + offset
                if name not in gene_regions:
                    gene_regions[name] = entry
                elif gene_regions[name][2] - gene_regions[name][1] < entry[2] - entry[1]:
                    gene_regions[name] = entry
    return gene_regions, unknown_genes


def build_seq_dict_names(gene_regions, assembly):
    """
    Takes gene_regions dict from names_to_seqs and builds a new dict with actual sequences
    uses call_twoBitToFa - returns a dict in the form {<gene>:<sequence>}
    """
    gene_sequences = dict()
    for gene in gene_regions:
        chrom, start, end, name, exon_count, strand = gene_regions[gene]
        seq = call_twoBitToFa(chrom, start, end, assembly)
        gene_sequences[gene] = seq
    return gene_sequences


def callHgsql(database, command):
    """
    Run hgsql command using subprocess, return stdout data if no error.
    Code from genecats (ENCODE)
    """
    cmd = ["hgsql", database, "-Ne", command]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmdout, cmderr = p.communicate()
    if p.returncode != 0:
        # keep command arguments nicely quoted
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return cmdout


def call_twoBitToFa(chrom, start, end, assembly):
    """
    Calls twoBitToFa on specific regions. Returns concatenated sequence
    """
    genome_path = "/cluster/data/{}/{}.2bit".format(assembly, assembly)
    cmd = ["twoBitToFa", "-seq=" + chrom, "-start=" + str(start), "-end=" + str(end), genome_path, "/dev/stdout"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmdout, cmderr = p.communicate()
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return "".join(cmdout.split()[1:])

def main(args):
    args = parse_args(args)
    #use hgsql to find gene regions of genes of interest - store failed sql calls
    gene_regions, unknown_genes = names_to_seqs(args.names, args.assembly, args.offset)
    #parse sql results into a dict (and remove duplicates retaining longest hit)
    gene_sequences = build_seq_dict_names(gene_regions, args.assembly)
    for gene in gene_sequences:
        args.fasta.write(">{}\n{}\n".format(gene, gene_sequences[gene]))
    for gene in unknown_genes:
        args.failed.write("{}\n".format(gene))


if __name__ == "__main__":
    sys.exit(main(sys.argv))


