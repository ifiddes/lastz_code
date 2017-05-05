
# coding: utf-8

# In[15]:

import pysam
import os
import pandas as pd
from pyfaidx import Fasta
from collections import *
from tools.intervals import *
from tools.misc import *
from tools.procOps import *
from tools.fileOps import *
from tools.bio import *
from tools.psl import *
from itertools import *
import bisect


# In[5]:

# first, construct a map of sequence positions to alignment positions
aln_f = Fasta('notch2nl_alignment.fa')
seq_aln_map = defaultdict(dict)
for name, seq in aln_f.iteritems():
    seq_pos = 0
    for aln_pos, x in enumerate(str(seq)):
        seq_aln_map[name][seq_pos] = aln_pos
        if x != '-':
            seq_pos += 1


# In[182]:

# find maximum position for reversing negative strand
max_pos = {x: max(y.keys()) for x, y in seq_aln_map.iteritems()}


# In[193]:

# next, construct a map of hg38 positions to sequence positions using the alignment
hg38_map = {}
for rec in pysam.Samfile('hg38_mapped.bam'):
    m = {y: x for x, y in rec.aligned_pairs}
    if rec.qname in ['NOTCH2', 'NOTCH2NL-A', 'NOTCH2NL-B']:
        m = {x: max_pos[rec.qname] - y for x, y in m.iteritems()}
    hg38_map[rec.qname] = m


# In[194]:

# construct a table mapping each alignment position to all hg38 positions
r = defaultdict(dict)
for name, pos_map in hg38_map.iteritems():
    for hg38_pos, seq_pos in pos_map.iteritems():
        aln_pos = seq_aln_map[name][seq_pos]
        r[name][aln_pos] = hg38_pos

df = pd.DataFrame.from_dict(r)
df.head()


# In[196]:

# now invert this map, so that we have our hg38 -> aln map
final_map = {}
for name in r:
    for aln_pos in r[name]:
        hg38_pos = r[name][aln_pos]
        assert hg38_pos not in final_map
        final_map[hg38_pos] = aln_pos


# In[197]:

# load the intervals we are interested in
start_stop_positions = {x.split()[3]: (x.split()[1], x.split()[2]) for x in open('n2_regions.bed')}
start_stop_positions = {x: map(int, y) for x, y in start_stop_positions.iteritems()}

# construct our regions of interest
regions_of_interest = []
for start, stop in start_stop_positions.itervalues():
    regions_of_interest.append(ChromosomeInterval('chr1', start, stop, '.'))


# In[198]:

# load all the alignments, by qname
aln_map = {}
base_dir = '/hive/users/ifiddes/notch2nl_berkeley_data/E2del19N_E2del68_combined_longranger/E2del68_E2del19N_combined/new-assembly/hg38_mapped_contigs'
for name in ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2', 'N1', 'N2']:
    s = os.path.join(base_dir, '{}.hg38.bam'.format(name))
    aln_map[name] = list(pysam.Samfile(s))


# In[199]:

# extract all alignments, discarding those not in our ROIs
blocks = {}
for name in aln_map:
    b = []
    for aln in aln_map[name]:
        if aln.is_unmapped:
            continue
        start = aln.reference_start
        end = aln.reference_end
        c = ChromosomeInterval('chr1', start, end, '.')
        if not interval_not_intersect_intervals(regions_of_interest, c) and len(c) > 100:
            b.append([start, end, aln.qname, aln.seq, aln.qstart, aln.qstart + aln.alen])
    blocks[name] = b


# In[200]:

for b in blocks['A1']:
    print b[0], b[1], b[2], b[4], b[5]


# In[249]:

# for each alignment, determine the interval in the MSA
# use a closest finding algorithm to deal with unaligned regions

class Record(object):
    def __init__(self, msa_start, msa_stop, qname, seq, seq_start, seq_end):
        self.msa_start = msa_start
        self.msa_stop = msa_stop
        self.qname = qname
        self.seq = seq
        self.seq_start = seq_start
        self.seq_end = seq_end

def find_closest(numeric_list, query_number):
    """
    Given a list of numbers, and a single query number, find the number in the sorted list that is numerically
    closest to the query number. Uses list bisection to do so, and so should be O(log n)
    """
    sorted_numeric_list = sorted(numeric_list)
    pos = bisect.bisect_left(sorted_numeric_list, query_number)
    if pos == 0:
        return sorted_numeric_list[0]
    if pos == len(sorted_numeric_list):
        return sorted_numeric_list[-1]
    before = sorted_numeric_list[pos - 1]
    after = sorted_numeric_list[pos]
    if after - query_number < query_number - before:
        return after
    else:
        return before

sorted_positions = sorted(final_map.keys())

msa_blocks = {}
for name, b in blocks.iteritems():
    mb = []
    for start, end, qname, seq, seq_start, seq_end in b:
        closest_start = find_closest(sorted_positions, start)
        closest_stop = find_closest(sorted_positions, end)
        msa_start = final_map[closest_start]
        msa_stop = final_map[closest_stop]
        if msa_start > msa_stop:  # handle negative strand
            msa_start, msa_stop = msa_stop, msa_start
        mb.append(Record(msa_start, msa_stop, qname, seq, seq_start, seq_end))
    # sort these to be in notch2nl order
    mb = sorted(mb, key=lambda x: x.msa_start)
    msa_blocks[name] = mb


# In[267]:

for x in msa_blocks['A2']:
    print x.msa_start, x.msa_stop, x.qname, x.seq_start, x.seq_end, len(x.seq)


# In[261]:

# filter blocks for those that are entirely a subset of another
# this removes misassemblies
filtered_blocks = {}
for name, mb in msa_blocks.iteritems():
    intervals = [ChromosomeInterval('', x.msa_start, x.msa_stop, '.', x) for x in mb]
    bad_intervals = set()
    for i1, i2 in combinations(intervals, 2):
        if i1.proper_subset(i2):
            bad_intervals.add(i1)
        elif i2.proper_subset(i1):
            bad_intervals.add(i2)
    filtered_blocks[name] = [x.data for x in intervals if x not in bad_intervals]


# In[266]:

for x in filtered_blocks['A2']:
    print x.msa_start, x.msa_stop, x.qname, x.seq_start, x.seq_end, len(x.seq)


# In[251]:

def find_overlap(b1, b2, spacer=50):
    delta = b1.msa_stop - b2.msa_start 
    b1_seq = b1.seq[-(delta + spacer):]
    b2_seq = b2.seq[:delta + spacer]
    aln = perform_aln(b1_seq, b2_seq)
    if len(aln) == 0:
        s = b1.seq[:-delta]
        assert len(s) > 0
        return s
    else:
        psl = sorted([PslRow(x.split('\t')) for x in aln], key=lambda x: x.coverage)[-1]
        cutoff = psl.query_coordinate_to_target(psl.q_end - 1) + 1
        return b1.seq[:-cutoff]

def merge_same(b1, b2):
    with TemporaryFilePath() as f, TemporaryFilePath() as f2:
        with open(f, 'w') as outf:
            write_fasta(outf, 'b1', b1.seq)
            write_fasta(outf, 'b2', b2.seq)
        cmd = ['muscle', '-in', f, '-out', f2]
        run_proc(cmd)
        fa = Fasta(f2)
        seq = []
        for x, y in zip(*[fa['b1'], fa['b2']]):
            if x == '-':
                seq.append(y)
            else:
                seq.append(x)
        return ''.join(seq)
    
def perform_aln(b1_seq, b2_seq):
    with TemporaryFilePath() as b1_f, TemporaryFilePath() as b2_f:
        with open(b1_f, 'w') as b1_f_h:
            write_fasta(b1_f_h, 'b1', b1_seq)
        with open(b2_f, 'w') as b2_f_h:
            write_fasta(b2_f_h, 'b1', b2_seq)
        cmd = ['blat', b1_f, b2_f, '-noHead', '/dev/stdout']
        return call_proc_lines(cmd)[:-2]


# In[278]:

seqs = {}
for name, mb in filtered_blocks.iteritems():
    seq = []
    a, b = tee(mb)
    _ = next(b, None)
    merged = set()
    for i, (b1, b2) in enumerate(izip(a, b)):
        if b1.qname in merged:
            continue
        #print i, b1.qname, b2.qname
        if b1.qname == b2.qname:
            #print 'merge', len(m)
            m = merge_same(b1, b2)
            seq.append(m)
            merged.add(b1.qname)
        elif b2.msa_start < b1.msa_stop:  # we have an overlap, try to resolve via pairwise alignment
            o = find_overlap(b1, b2)
            #print 'overlap', len(o)
            seq.append(o)
        elif b1.msa_stop == b2.msa_start:
            #print 'same start/stop', len(b1.seq)
            seq.append(b1.seq)
        else:
            seq.append(b1.seq)
            seq.append(''.join(['N'] * 100))
            #print 'gap', len(b1.seq), delta
    seq.append(b2.seq)  # add the last sequence
    seqs[name] = seq


# In[279]:

with open('a1.contigs.fa', 'w') as outf:
    for i, x in enumerate(seqs['A1']):
        write_fasta(outf, str(i), x)


# In[280]:

get_ipython().magic(u'mkdir -p assemblies_v2')
from tools.procOps import *
for name, seq in seqs.iteritems():
    with open('assemblies_v2/{}.fa'.format(name), 'w') as outf:
        write_fasta(outf, name, ''.join(seq))
    cmd = [['bwa', 'mem', '-t', '4', '-x', 'intractg', '/hive/groups/recon/10x_genomics/references/refdata-hsapiens-hg38/fasta/genome.fa', 
            'assemblies_v2/{}.fa'.format(name)],
           ['samtools', 'view', '-b', '-'],
           ['sambamba', 'sort', '-o', 'assemblies_v2/{}.hg38.bam'.format(name), '/dev/stdin']]
    run_proc(cmd)


# In[ ]:



