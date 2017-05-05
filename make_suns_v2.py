#!/usr/bin/env python

"""
suns to vcf (hg19)
"""
import sys


rc = {"A":"T", "T":"A", "G":"C", "C":"G"}

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



aligned = {x[0].replace("_extraction","").replace("_(reversed)",""): x[2] for x in read_fasta(open(sys.argv[1]))}

#offset = 120551343 #for use with good_region_realign_hg19.fasta
offset = 120535415 #for use with 


suns = {"A":[],"B":[],"C":[],"D":[],"AB":[],"N":[]}
gapN = 0
gapA, gapB, gapC, gapD = 0, 0, 0, 0
#ref base, alt base, hg19 pos, alignment pos, seq_pos
for i in xrange(len(aligned["notch2"])):
    n,b,a,d,c = aligned['notch2'][i], aligned['notch2nl-B'][i], aligned['notch2nl-A'][i], aligned['notch2nl-D'][i], aligned['notch2nl-C'][i]
    if n != "-" and b != "-" and a != "-" and d != "-" and c != "-":
        if n != a and n != b and n != d and n != c:
            suns["N"].append([n, a, i-gapN+offset, i, i-gapN]) #n2 has whatever base the others have
        elif a != n and a != b and a != c and a != d:
            suns["A"].append([n, a, i-gapN+offset, i, i-gapA])
        elif b != n and b != a and b != c and b != d:
            suns["B"].append([n, b, i-gapN+offset, i, i-gapB])
        elif c != n and c != a and c != b and c != d:
            suns["C"].append([n, c, i-gapN+offset, i, i-gapC])
        elif d != n and d != c and d != b and d != b:
            suns["D"].append([n, d, i-gapN+offset, i, i-gapD])
        elif a == b and a != d and a != c and a != n:
            suns["AB"].append([n, a, i-gapN+offset, i, i])
    if n == "-":
        gapN += 1
    if a == "-":
        gapA += 1
    if b == "-":
        gapB += 1
    if c == "-":
        gapC += 1
    if d == "-":
        gapD += 1

header = "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


QUAL="200"; FILTER="."; INFO="."; ID="."; chrom="chr1"

for para in suns:
    outf = open(para + ".vcf", "w")
    outf.write(header)
    for ref, alt, pos, aln_pos, seq_pos in suns[para]:
        line = "\t".join(map(str,[chrom, pos, ID, ref, alt, QUAL, FILTER, INFO])) + "\n"
        outf.write(line)
    outf.close()

wl = open("whitelist.txt","w")
all_positions = []
for para in suns:
    for ref, alt, pos, aln_pos, seq_pos in suns[para]:
        all_positions.append([para, pos])

for x in all_positions:
    wl.write("\t".join(map(str,x))); wl.write("\n")

wl.close()

#make map to CHM1
#starts = {"A":(147536494-243+15-1,"-"),"B":(149987434-228-1,"-"),"C":(150820087+230-1,"+"),"D":(121442729+236-1,"+"),"N":(120666926-112,"-")}
starts = {"A":(147520393, "-"), "B":(149971335, "-"), "C":(150786389, "+"), "D":(121413083, "+"), "N":(120650883, "-")}
chm1 = []
for para in starts:
    for ref, alt, pos, aln_pos, seq_pos in suns[para]:
        offset, strand = starts[para]
        if strand == "-":
            new_pos = seq_pos + offset
        else:
            new_pos = offset - seq_pos
        chm1.append(["chr1",new_pos-1,new_pos,para + "_"+str(pos)])

outf = open("n2nl_map_CHM1.bed","w")
outf.write("track name=CHM1_map\n")
for x in chm1:
    outf.write("\t".join(map(str,x))); outf.write("\n")