import pysam, os, sys, argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", "-r", type=argparse.FileType("r"), 
        help="Reference fasta file")
    parser.add_argument("--sam", "-s", type=str, help="Samfile")
    parser.add_argument("--out", "-o", type=argparse.FileType("w"),
        help="output MAF-style-ish file.")
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


def convert_sam_record(record, target):
    """Takes a single record from a sam and a str representing a sequence
    and converts this toper a aligned fasta (MAF) style format
    """
    read, ref = list(), list()
    pos, ref_pos = 0, record.pos
    for oper, num in record.cigar:
        if oper == 0 or oper == 7:
            #treating M and = equally
            read.append(record.seq[pos : pos + num])
            ref.append(record.seq[pos : pos + num])
            pos += num; ref_pos += num
        elif oper == 1:
            #insertion
            read.append(record.seq[pos : pos + num])
            ref.append("-" * num)
            pos += num
        elif oper == 2 or oper == 3:
            #treating D and N equally
            read.append("-" * num)
            ref.append(target[ref_pos : ref_pos + num])
            ref_pos += num
        elif oper == 8:
            #mismatch
            read.append(record.seq[pos : pos + num])
            ref.append(target[ref_pos : ref_pos + num])
            pos += num; ref_pos += num
        elif oper == 4:
            #soft clip means ignore that part of read
            pos += num
        elif oper == 5:
            #hard-clipped sequences can be ignored here
            continue
        elif oper == 6:
            #ignore padding as well
            continue
    return ("".join(read).upper(), "".join(ref).upper())

def main(args):
    args = parse_args(args)
    referenceDict = {name:seq for name, com, seq in read_fasta(args.reference)}
    sam = pysam.Samfile(args.sam, "r")
    for record in sam:
        rname = sam.getrname(record.tid)
        qname = record.qname
        rseq = referenceDict[rname]
        seq, ref = convert_sam_record(record, rseq)
        args.out.write(">{}\n{}\n>{}\n{}\n".format(qname,seq,rname,ref))

if __name__ == '__main__':
    sys.exit(main(sys.argv))
