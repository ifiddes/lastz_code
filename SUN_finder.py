#!/usr/bin/env python

from Bio import SeqIO
import sys, os, subprocess
#expects that geneious was used to reference assemble to N2
#expects the first argument to be a file listing all of these fastas
with open(sys.argv[1]) as f:
    SUNs = {}
    for fasta in f:
        records = [x for x in SeqIO.parse(fasta.rstrip(),"fasta")]
        for r in records:
            if r.id not in SUNs:
                SUNs[r.id] = []
            if r.id == "Notch2":
                offset = int(r.description.split(" ")[-1])
                tmp = r

        n2_pos = -1
        for aln_pos, n2_base in enumerate(str(r.seq)):
            #iterate over the position and base in Notch2
            if n2_base != "-":
                #n2_pos stores the position in the n2 sequence
                n2_pos += 1
            #list of bases at this alignment position
            bases = [str(x.seq[aln_pos]) for x in records]
            #definition of a SUN
            if "-" not in bases and len(set(bases)) == 2:
                positions = []
                for base in set(bases):
                    positions.append([a for a, x in enumerate(bases) if x == base])
                for x in xrange(len(positions)):
                    if len(positions[x]) == 1:
                        try:
                            alt = bases[positions[x-1][0]]
                        except:
                            alt = bases[positions[x+1][0]]
                        SUNs[records[positions[x][0]].id].append([n2_pos+offset, bases[positions[x][0]], alt, len(records)])               

for s in SUNs:
    SUNs[s] = sorted(SUNs[s], key = lambda x:x[0])


with open("whitelist.txt", "w") as outf:
    outf.write("#paralog\thg38_n2_position\tSUN\talt\texp_frac\n")
    for para, suns in SUNs.iteritems():
        for sun in suns:
            outf.write("\t".join([para,"\t".join(map(str,sun))])+"\n")

for para, suns in SUNs.iteritems():
    if para == "Notch2":
        para = "N"
    else:
        para = para.split("-")[-1]
    with open("{}.vcf".format(para),"w") as outf:
        outf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
        for pos, ref, alt, frac in suns:
            outf.write("\t".join(map(str,["chr1",pos-1,".",ref,alt,"200",".","."]))+"\n")
    subprocess.call(["bgzip", "{}.vcf".format(para)])
    subprocess.call(["tabix", "-p", "vcf", "{}.vcf.gz".format(para)])
    with open("{}.bedGraph".format(para), "w") as outf:
        outf.write("track name={}\n".format(para))
        for pos, ref, alt, frac in suns:
            outf.write("\t".join(["chr1",str(pos-1),str(pos),alt])+"\n")