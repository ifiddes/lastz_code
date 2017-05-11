import argparse
from tools.procOps import *
from tools.intervals import *
from tools.fileOps import *
from tools.bio import *
import pysam
from collections import defaultdict
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--a-fastas', nargs=2, required=True)
    parser.add_argument('--a-ref', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-A.fasta')
    parser.add_argument('--a-bed', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-A.bed')
    parser.add_argument('--b-fastas', nargs=2, required=True)
    parser.add_argument('--b-ref', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-B.fasta')
    parser.add_argument('--b-bed', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-B.bed')
    parser.add_argument('--c-fastas', nargs=2, required=True)
    parser.add_argument('--c-ref', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-C.fasta')
    parser.add_argument('--c-bed', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-C.bed')
    parser.add_argument('--d-fastas', nargs=2, required=True)
    parser.add_argument('--d-ref', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-D.fasta')
    parser.add_argument('--d-bed', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2NL-D.bed')
    parser.add_argument('--n-fastas', nargs=2, required=True)
    parser.add_argument('--n-ref', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2.fasta')
    parser.add_argument('--n-bed', default='/hive/users/ifiddes/notch2nl_berkeley_data/hg38_haplotype_indices/Notch2.bed')
    parser.add_argument('--out-fasta', required=True)
    return parser.parse_args()


def map_to_reference(fa, reference):
    bam_path = get_tmp_file()
    cmd = [['bwa', 'mem', '-t', '4', '-x', 'intractg', reference, fa],
           ['samtools', 'view', '-b', '-'],
           ['sambamba', 'sort', '-o', bam_path, '/dev/stdin']]
    run_proc(cmd)
    return bam_path


def parse_intervals(bed):
    return [ChromosomeInterval(x.split()[0], x.split()[1], x.split()[2], x.split()[-1], x.split()[3])
             for x in open(bed)]


def find_exons(bed, bam_path, intervals):
    results = defaultdict(list)
    for r in pysam.Samfile(bam_path):
        if not r.is_unmapped:
            for i in intervals:
                if r.reference_start <= i.start and r.reference_end >= i.stop:
                    read_positions, ref_positions = zip(*r.aligned_pairs)
                    start = read_positions[ref_positions.index(i.start)]
                    stop = read_positions[ref_positions.index(i.stop)]
                    results[i.data].append(r.seq[start:stop])
    #assert len(results) == 5, (intervals, results.keys(), bam_path)
    if len(results) != 5:
        missing_exons = {'E1', 'E2', 'E3', 'E4', 'E5'} - set(results.viewkeys())
        'Missing exons {}'.format(','.join(missing_exons))
    for name, r in results.iteritems():
        if len(r) > 1:
            print 'Found {} alignments for {}-{}.'.format(len(r), bam_path, name)
            if len(set(r)) == 1:
                print 'Sequences match'
            else:
                assert False

    return ''.join([results['E1'][0], results['E2'][0], results['E3'][0], results['E4'][0], results['E5'][0]])


def main():
    args = parse_args()
    with open(args.out_fasta, 'w') as outf:
        for bin_group, ref, bed, name in [[args.a_fastas, args.a_ref, args.a_bed, 'NOTCH2NL-A'],
                                          [args.b_fastas, args.b_ref, args.b_bed, 'NOTCH2NL-B'],
                                          [args.c_fastas, args.c_ref, args.c_bed, 'NOTCH2NL-C'],
                                          [args.d_fastas, args.d_ref, args.d_bed, 'NOTCH2NL-D'],
                                          [args.n_fastas, args.n_ref, args.n_bed, 'NOTCH2']]:
            intervals = parse_intervals(bed)
            for i, fa in enumerate(bin_group, 1):
                bam_path = map_to_reference(fa, ref)
                seq = find_exons(bed, bam_path, intervals)
                write_fasta(outf, '{}-{}'.format(name, i), seq)
                os.remove(bam_path)


if __name__ == '__main__':
    main()