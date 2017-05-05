outA = open("n2nl_chm1_A.vcf", "w")
outB = open("n2nl_chm1_B.vcf", "w")
outC = open("n2nl_chm1_C.vcf", "w")
outD = open("n2nl_chm1_D.vcf", "w")
outN = open("n2nl_chm1_N.vcf", "w")
outAB = open("n2nl_chm1_AB.vcf", "w")

for line in open(sys.argv[1]):
	line = line.split()
	paralog = line[3].split("_")[0]
	if paralog == "A":
		outA.write("\t".join(line)); outA.write("\n")
	elif paralog == "B":
		outB.write("\t".join(line)); outB.write("\n")
	elif paralog = "C":
		outC.write("\t".join(line)); outC.write("\n")
	elif paralog == "D":
		outD.write("\t".join(line)); outD.write("\n")
	elif paralog == "N":
		outN.write("\t".join(line)); outN.write("\n")
	elif paralog == "AB(A)":
		line[3] = line[3].replace("(A)", "")
		outAB.write("\t".join(line)); outAb.write("\n")

outA.close(); outB.close(); outC.close(); outD.close(); outN.close(); outAB.close()