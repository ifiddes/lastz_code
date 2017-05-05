"""
Simple script to find the coverage of a BED given a genome
"""
import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('bed')
    parser.add_argument('chromSizes')
    return parser.parse_args()


def format_ratio(numerator, denominator):
    """
    Convenience function that converts two numbers, integer or no, to a ratio
    """
    if denominator == 0:
        return float("nan")
    return float(numerator) / denominator


def parse_bed(bed):
    c = 0
    for l in open(bed):
        l = l.split()
        c += int(l[2]) - int(l[1])
    return c


def parse_sizes(sizes):
    c = 0
    for l in open(sizes):
        l = l.split()
        c += int(l[1])
    return c


def main():
    args = parse_args()
    c = parse_bed(args.bed)
    s = parse_sizes(args.chromSizes)
    print 'Percent coverage: {:.3%}'.format(format_ratio(c, s))


if __name__ == '__main__':
    main()