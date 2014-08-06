#!/bin/bash

GENOME=`basename $1 | cut -d "." -f 1`
TARGET="/hive/data/genomes/$GENOME/$GENOME.2bit"
QUERY="/cluster/home/ifiddes/all_ZNF_FIXED.2bit"

TSIZES="/hive/data/genomes/$GENOME/chrom.sizes"
QSIZES="/cluster/home/ifiddes/all_ZNF_FIXED.sizes"

mkdir $GENOME

#parse axt and chain
python /cluster/home/ifiddes/code/fix_axt_align_seq2_thing.py $1 | \
axtChain -verbose=0 -scoreScheme=/scratch/data/blastz/human_chimp.v2.q  \
-minScore=5000 -linearGap=medium stdin $TARGET $QUERY stdout | \
chainAntiRepeat $TARGET $QUERY stdin $GENOME/$GENOME.chains

#make nets
chainPreNet $GENOME/$GENOME.chains $TSIZES $QSIZES stdout | \
chainNet stdin -minSpace=1 $TSIZES $QSIZES stdout /dev/null | \
netSyntenic stdin $GENOME/$GENOME.noClass.net

#generate liftover chain
netChainSubset -verbose=0 $GENOME/$GENOME.noClass.net $GENOME/$GENOME.chains \
stdout | chainStitchId stdin $GENOME/$GENOME.liftover.chain

#swap ref to ZNF
chainStitchId $GENOME/$GENOME.liftover.chain stdout | chainSwap stdin \
stdout | chainSort stdin $GENOME/$GENOME.tBest.chain

#net these to get reciprocal best net
chainPreNet $GENOME/$GENOME.tBest.chain $QSIZES $TSIZES stdout | \
chainNet -minSpace=1 -minScore=0 stdin $QSIZES $TSIZES stdout /dev/null | \
netSyntenic stdin $GENOME/$GENOME.ZNFref.rBest.net

#extract reciprocal best chain; ZNF ref'd
netChainSubset $GENOME/$GENOME.ZNFref.rBest.net $GENOME/$GENOME.tBest.chain \
stdout | chainStitchId stdin $GENOME/$GENOME.ZNFref.rBest.chain

#swap to get $GENOME ref'd best chain
chainSwap $GENOME/$GENOME.ZNFref.rBest.chain stdout | \
chainSort stdin $GENOME/$GENOME.rBest.chain

#net to get $GENOME ref'd reciprocal best net
chainPreNet $GENOME/$GENOME.rBest.chain $TSIZES $QSIZES stdout | \
chainNet -minSpace=1 -minScore=0 stdin $TSIZES $QSIZES stdout /dev/null | \
netSyntenic stdin $GENOME/$GENOME.rBest.net

#generate bed files of both nets
netToBed -maxGap=1 $GENOME/$GENOME.rBest.net $GENOME/$GENOME.rBest.bed
netToBed -maxGap=1 $GENOME/$GENOME.ZNFref.rBest.net \
$GENOME/$GENOME.ZNFref.rBest.bed

#generate psls from chains
chainToPsl $GENOME/$GENOME.rBest.chain $TSIZES $QSIZES $TARGET $QUERY \
$GENOME/$GENOME.rBest.psl
chainToPsl $GENOME/$GENOME.ZNFref.rBest.chain $QSIZES $TSIZES $QUERY $TARGET \
$GENOME/$GENOME.ZNFref.rBest.psl