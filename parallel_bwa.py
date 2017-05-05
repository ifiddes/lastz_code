import sys
import argparse
import os
os.environ['PYTHONPATH'] = './:/hive/users/ifiddes/ihategit/pipeline/:/hive/users/ifiddes/ihategit/pipeline/submodules:/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio'
sys.path.extend(['./', '/hive/users/ifiddes/ihategit/pipeline/', '/hive/users/ifiddes/ihategit/pipeline/submodules', '/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio', '/hive/users/ifiddes/comparativeAnnotator'])
from jobTree.scriptTree.stack import Stack, Target
from pycbio.sys.fileOps import tmpFileGet
from pycbio.sys.procOps import runProc
from pycbio.sys.dataOps import grouper
from sonLib.bioio import fastaRead, TempFileTree


def run_bwa(target, ref_idx, tmp_fasta, out_dir, num_threads):
    tmp_sort = tmpFileGet()
    out_path = tmpFileGet(tmpDir=out_dir)
    cmd = [['bwa', 'mem', '-t', num_threads, ref_idx, tmp_fasta],
           ['samtools', 'view', '-b', '-'],
           ['samtools', 'sort', '-O', 'bam', '-T', tmp_sort, '-']]
    runProc(cmd, stdout=out_path)


def run_bwa_wrapper(target, args):
    for base_path, dirs, files in os.walk(args.splitFastaDir):
        if files:
            for f in files:
                p = os.path.join(base_path, f)
                target.addChildTargetFn(run_bwa, memory=8 * 1024 ** 3,
                                        args=(args.reference, p, target.getGlobalTempDir(), args.defaultCpu))
    target.setFollowOnTargetFn(cat, args=(args,))


def cat(target, args):
    fofn = tmpFileGet()
    files = [os.path.join(target.getGlobalTempDir(), x) for x in os.listdir(target.getGlobalTempDir())]
    files = [x for x in files if os.path.isfile(x)]
    assert len(files) > 0
    with open(fofn, 'w') as outf:
        for x in files:
            outf.write(x + "\n")
    cmd = ['samtools', 'merge', '-b', fofn, args.outBam]
    runProc(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outBam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--splitFastaDir", required=True)
    parser.add_argument('--numThreads', default=4)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    args.defaultCpu = args.numThreads
    i = Stack(Target.makeTargetFn(run_bwa_wrapper, memory=8 * 1024 ** 3, args=(args,))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from parallel_bwa import *
    main()
