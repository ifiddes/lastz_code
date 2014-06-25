#!/usr/bin/bash

while read g; do
    python2.7 ~/code/lastz_32_general_format.py --target /cluster/data/$g/$g.2bit --cores 10 --chain --fasta ~/all_ZNF_FIXED.fasta > $g.chains
    python2.7 ~/code/merge_chains.py $g.chains > $g.chains.merged
    python2.7 ~/code/merged_to_bed.py $g.chains.merged > $g.chains.merged.bed
done < $1