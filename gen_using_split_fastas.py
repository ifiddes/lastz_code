import sys, os, itertools, argparse, subprocess, pipes, glob
from multiprocessing import Pool

def parse_args(args):
    """parses arguments using argparse"""
    parser = argparse.ArgumentParser(prog="pipeline",description=__doc__, usage= \
            "pipeline [args]")
    parser.add_argument("--cores", default=1, type=int, help="Number of cores (default: 1)")
    parser.add_argument("--fastdir", required=True, help="Genome fasta directory")
    parser.add_argument("--znfdir", required=True, help="ZNF fasta directory")
    args = parser.parse_args()
    return args


def call_lastz_32_chain(call):
    """
    Generic function that runs lastz with --chain on a string
    And a path pointing to a 1-entry fasta file containing a query chr_seq
    """
    z, g = call
    cmd = ["lastz_32", "--format=general:name1,zstart1,end1,name2,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx-", "--chain", "--ambiguous=n", z, g]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    cmdout, cmderr = p.communicate()
    if p.returncode != 0:
        cmdstr = " ".join([pipes.quote(arg) for arg in cmd])
        raise Exception("Error from: " + cmdstr + ": " + cmderr)
    return cmdout


def main(args):
	args = parse_args(args)
	z = glob.glob(args.znfdir + "/*")
	g = glob.glob(args.fastdir + "/*")
	pool = Pool(processes = args.cores)
	i = itertools.product(z, g)
	r = pool.map(call_lastz_32_chain, i)
	pool.close()
	pool.join()
	for x in r:
		sys.stdout.write(x)



if __name__ == "__main__":
    sys.exit(main(sys.argv))
    




