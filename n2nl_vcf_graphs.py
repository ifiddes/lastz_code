#!/usr/bin/env python
"""
Parse vcfs from pipeline and create plots across the five genes
"""
import vcf, sys, os, argparse, subprocess, time
from itertools import izip
from collections import OrderedDict as od

def parse_args(args):
    parser = argparse.ArgumentParser(prog="graphVcf",description=__doc__)
    parser.add_argument("--folder", type=str, required=True, help="folder containing samples in form <ID>/<A,B,C,D>.vcf")
    parser.add_argument("--tsv", type=str, required=True, help="file to write full tsv to")
    parser.add_argument("--summary", type=str, required=True, help="file to summary tsv to (median/average/variance)")
    parser.add_argument("--graph", type=str, required=True, help="root of file to write plots")
    parser.add_argument("--whitelist", type=argparse.FileType("r"), required=True, help="whitelist txt file (tsv, 1st field = paralog, 2nd field = genome pos")
    return parser.parse_args()


def make_R(vcf_file, whitelist, start, stop):
    """Function that returns a list with every position represented (for barplot)"""
    vcf_reader = vcf.Reader(file(vcf_file))
    results = dict()
    for record in vcf_reader:
        if record.POS in whitelist:
            results[record.POS] = float(record.INFO["ALTFRAC"][0])
    final = list()
    for i in xrange(start, stop):
        if i in results:
            final.append(results[i])
        else:
            final.append(0)
    return final


def parse_whitelist(whitelist):
    """Parses a whitelist txt file into a dict mapping a paralog to a list of positions"""
    wl = dict()
    for line in whitelist:
        line = line.split()
        paralog, pos = line[:2]
        if paralog not in wl:
            wl[paralog] = list()
        wl[paralog].append(int(pos))
    return wl

def mean(ll):
    """per column mean of a list of lists"""
    return [sum(l) / len(l) * 1.0 for l in map(list, zip(*ll))]

def median(ll):
    """per column median of a list of lists"""
    s = [sorted(l) for l in map(list, zip(*ll))]
    return [(x[len(x) / 2] + x[len(x) / 2 - 1]) / 2.0 if not (len(x) % 2) else x[len(x) / 2] for x in s]

def var(ll):
    """per column variance of a list of lists"""
    ll = map(list, zip(*ll))
    return [sum((a - val) ** 2 for val in l) / len(l) for a, l in izip([sum(x) / len(x) for x in ll], ll)]
        

def main(args):
    paralogs = ["AB", "A", "B", "C", "D", "N"]
    start = 120535758
    stop = 120598753

    args = parse_args(args)
    wl = parse_whitelist(args.whitelist)

    outf = open(args.tsv, "w")

    samples = sorted(os.listdir(args.folder))
    sample_paths = [os.path.join(args.folder, x) for x in samples]
    names = [os.path.basename(x) for x in samples]

    final = od()

    for sample, n in izip(sample_paths, names):
        combined = od()
        vcfs = [os.path.join(sample, x) for x in ["AB.vcf", "A.vcf", "B.vcf", "C.vcf", "D.vcf", "N.vcf"]]
        for v, p in izip(vcfs, paralogs):
            combined[p] = make_R(v, wl[p], start, stop)
            final[n] = combined

    for sample, paralog_dict in final.iteritems():
        for paralog, csv in paralog_dict.iteritems():
            outf.write("\t".join(map(str, [sample + "_" + paralog] + csv)))
            outf.write("\n")

    outf.close()

    call = ["Rscript", "/cluster/home/ifiddes/code/n2nl_barplots.R", args.tsv, args.graph + "_full.png"]
    p1 = subprocess.Popen(call)

    outf = open(args.summary, "w")
    rearranged=od()

    for sample, paralog_dict in final.iteritems():
        for paralog, csv in paralog_dict.iteritems():
            if paralog not in rearranged:
                rearranged[paralog] = list()
            rearranged[paralog].append(csv)

    for paralog, csvs in rearranged.iteritems():
        outf.write("\t".join(map(str, ["median", paralog] + median(csvs)))); outf.write("\n")
        outf.write("\t".join(map(str, ["mean", paralog] + mean(csvs)))); outf.write("\n")
        outf.write("\t".join(map(str, ["var", paralog] + var(csvs)))); outf.write("\n")

    outf.close()

    call = ["Rscript", "/cluster/home/ifiddes/code/n2nl_summary_barplots.R", args.summary, args.graph + "_summary.png"]
    p2 = subprocess.Popen(call)

    #while True:
    #    time.sleep(0.5)
    #    if p1.poll() and p1.poll():
    #        break


if __name__ == "__main__":
    sys.exit(main(sys.argv))
