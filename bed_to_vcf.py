#!/usr/bin/env python
"""script to convert a given bed file with associated fasta file to a simple vcf file
Currently hard coded on line 51 to handle the naming scheme in my notch2nl paralog analysis"""

import sys, os, argparse

header = "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"

def parse_args(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("--bed", type=argparse.FileType("r"), help="input bed file")
	parser.add_argument("--fasta", type=argparse.FileType("r"), help="input fasta file")
	parser.add_argument("--vcf", type=str, help="vcf files base")
	return parser.parse_args()

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


def read_bed(file_handle):
	for line in file_handle:
		line = line.rstrip().split("\t")
		if len(line) == 4:
			yield line
		else:
			continue


def main(args):
	args = parse_args(args)
	outA = open(args.vcf + "_A.vcf", "w")
	outB = open(args.vcf + "_B.vcf", "w")
	outC = open(args.vcf + "_C.vcf", "w")
	outD = open(args.vcf + "_D.vcf", "w")
	outN = open(args.vcf + "_N.vcf", "w")
	outAB = open(args.vcf + "_AB.vcf", "w")
	QUAL="200"; FILTER="."; INFO="."; ID="."
	fasta = dict([(x[0], x[2]) for x in read_fasta(args.fasta)])
	outA.write(header); outB.write(header); outC.write(header); outD.write(header)
	outN.write(header); outAB.write(header)
	for chrom, start, end, name in read_bed(args.bed):
		if chrom in fasta:
			pos = end #is this right?
			#paralog, hg19_pos, alt = name.split("_") #hg38 version
			paralog, alt = name.split("_") #hg19 version
			ref = fasta[chrom][int(end)].upper()
			line = "\t".join([chrom, pos, ID, ref, alt, QUAL, FILTER, INFO]) + "\n"
		if paralog == "A":
			outA.write(line)
		elif paralog == "B":
			outB.write(line)
		elif paralog == "C":
			outC.write(line)
		elif paralog == "D":
			outD.write(line)
		elif paralog == "N":
			outN.write(line)
		elif paralog == "AB(A)":
			outAB.write(line)

	outA.close(); outB.close(); outC.close(); outD.close(); outN.close(); outAB.close()		



if __name__ == "__main__":
    sys.exit(main(sys.argv))