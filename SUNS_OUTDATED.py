#!/usr/bin/env python

"""
iterate over aligned fasta from geneious and recreate SUN positions in region of interest.
Usage: <aligned fasta> <hg19 start pos (stranded!)
hg19 original pos: 120613523 CHM1 positions are currently hard coded
"""
import sys
aligned = [x.split() for x in open(sys.argv[1])]

rc = {"A":"T", "T":"A", "G":"C", "C":"G"}

seqs, names = [], []
for x in aligned:
	if x[0].startswith(">"):
		names.append(x[0].replace(">","").replace("_extraction",""))
	else:
		seqs.append(x[0])

#format: hg38 rel pos, hg19 rel pos, paralog
suns = {"A":[],"B":[],"C":[],"D":[],"AB":[],"N":[]}
gapN, gapA, gapB, gapC, gapD = 0, 0, 0, 0, 0
for i in xrange(len(seqs[0])):
	n,b,a,d,c = [x[i] for x in seqs]
	if n != "-" and b != "-" and a != "-" and d != "-" and c != "-":
		if n != a and n != b and n != d and n != c:
			suns["N"].append([i-gapN,i-gapN,rc[n]])
		elif a != n and a != b and a != c and a != d:
			suns["A"].append([i-gapA,i-gapN,rc[a]])
		elif b != n and b != a and b != c and b != d:
			suns["B"].append([i-gapB,i-gapN,rc[b]])
		elif c != n and c != a and c != b and c != d:
			suns["C"].append([i-gapC,i-gapN,c])
		elif d != n and d != c and d != b and d != b:
			suns["D"].append([i-gapD,i-gapN,d])
		elif a == b and a != d and a != c and a != n:
			suns["AB"].append([i-gapA,i-gapN,rc[a],i-gapB,i-gapN,rc[b]])
	if n == "-":
		gapN += 1
	elif a == "-":
		gapA += 1
	elif b == "-":
		gapB += 1
	elif c == "-":
		gapC += 1
	elif d == "-":
		gapD += 1

if len(sys.argv) == 2:
	hg19_ref_start = int(sys.argv[2])
else:
	hg19_ref_start = 120613523
#CHM1 = {"A":(147579521,147520745),"B":(150028031,149971687),"C":(150772632,150835803),"D":(121444749,121455782),"N":(120702621,120651235)}
#CHM1 = {"C":(150694367+65760, "+"), "B":(149971687, "-"), "A":(147520745, "-"), "D":(121444749, "+"), "N":(120651235,"-")}
final = []
CHM1 = {"C":(150772632-15000-645+472+2, "+"), "B":(149971687+77400+63, "-"), "A":(147520745+76800-25, "-"), "D":(121444749-64000+6, "+"), "N":(120651235+76931+840-8,"-")}


for paralog in suns:
	if paralog != "AB": #need to handle AB case separately
		for hg38_rel_pos, hg19_rel_pos, alt in suns[paralog]:
			hg19_pos = hg19_ref_start - hg19_rel_pos
			if CHM1[paralog][1] == "+":
				hg38_pos = CHM1[paralog][0] + hg38_rel_pos
			else: #paralog is on (-) strand
				hg38_pos = CHM1[paralog][0] - hg38_rel_pos
			if hg19_pos >= 120551500:
				final.append(map(str,["chr1",hg38_pos,hg38_pos+1,paralog+"_"+str(hg19_pos)+"_"+alt]))
	else:
		for a_pos, hg19_a_pos, a_alt, b_pos, hg19_b_pos, b_alt in suns[paralog]:
			hg19_pos = hg19_ref_start - hg19_a_pos
			chm1_a_pos = CHM1["A"][0] - a_pos
			chm1_b_pos = CHM1["B"][0] - b_pos
			final.append(map(str,["chr1",chm1_a_pos,chm1_a_pos+1,"AB(A)_"+str(hg19_pos)+"_"+a_alt]))
			final.append(map(str,["chr1",chm1_b_pos,chm1_b_pos+1,"AB(B)_"+str(hg19_pos)+"_"+b_alt]))

out_chm1 = open("SUNs_mapped_CHM1.bed","w")
out_chm1.write('track name="SUNS hg38"\n')
outf = open("SUNS_hg19_check.bed","w")
outf.write('track name="SUNs hg19 check"\n')
#write it in hg19 coords too for sanity tests
for line in final:
	out_chm1.write("\t".join(line)); out_chm1.write("\n")
	pos = line[-1].split("_")[1]
	line[-1] = line[-1].split("_")[0] + "_" + line[-1].split("_")[-1] #rewrote this for hg19 to save alt
	line[1] = line[2] = pos
	line[1] = str(int(line[2])-1)
	outf.write("\t".join(line));outf.write("\n")

outf.close()
out_chm1.close()



