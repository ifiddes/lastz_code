from __future__ import absolute_import
from argparse import ArgumentParser
import os
import logging
import random
import shutil
from tools.bio import *
from toil.job import Job
from toil.common import Toil
from tools.dataOps import grouper
from tools.fileOps import get_tmp_file
from tools.procOps import call_proc_lines

def setup(job, reference, target, chunk_size, ooc):
    ref_file_id = job.fileStore.writeGlobalFile(reference)
    tgt_file_id = job.fileStore.writeGlobalFile(target)
    if ooc is not None:
        ooc_file_id = job.fileStore.writeGlobalFile(ooc)
    else:
        ooc_file_id = None
    return job.addFollowOnJobFn(split, ref_file_id, tgt_file_id, chunk_size, ooc_file_id).rv()


def split(job, ref_file_id, tgt_file_id, chunk_size, ooc_file_id):
    rets = []
    tgt_file = job.fileStore.readGlobalFile(tgt_file_id)
    for chunk in grouper(read_fasta(tgt_file), chunk_size):
        rets.append(job.addChildJobFn(aln, chunk, ref_file_id, ooc_file_id, memory='8G').rv())
    return job.addFollowOnJobFn(combine, rets).rv()


def aln(job, chunk, ref_file_id, ooc_file_id):
    ref_file = job.fileStore.readGlobalFile(ref_file_id)
    if ooc_file_id is not None:
        ooc_file = job.fileStore.readGlobalFile(ooc_file_id)
    tmp = get_tmp_file()
    with open(tmp, 'w') as outf:
        for name, seq in chunk:
            write_fasta(outf, name, seq)
    cmd = ['blat', '-noHead', ref_file, tmp, '/dev/stdout']
    if ooc_file_id is not None:
        cmd.append('-ooc={}'.format(ooc_file))
    return call_proc_lines(cmd)[:-2]


def combine(job, rets):
    return '\n'.join(['\n'.join(x).rstrip() for x in rets]) + '\n'


def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--target', required=True)
    parser.add_argument('--chunk-size', default=500, type=int)
    parser.add_argument('--out-psl', required=True)
    parser.add_argument('--ooc')
    args = parser.parse_args()
    r = Job.Runner.startToil(Job.wrapJobFn(setup, os.path.abspath(args.reference),
                                           os.path.abspath(args.target), args.chunk_size, args.ooc,
                         memory='4G'), args)
    with open(args.out_psl, 'w') as outf:
        outf.write(r)


if __name__ == '__main__':
    main()