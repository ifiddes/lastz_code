#!/bin/bash

GENOME=`basename $1 | cut -d "." -f 1`
TARGET="/hive/data/genomes/$GENOME/$GENOME.2bit"
QUERY="/cluster/home/ifiddes/all_ZNF_FIXED.2bit"
TARGET_FASTA="/cluster/home/ifiddes/fastas/$GENOME.fa"
QUERY_FASTA="/cluster/home/ifiddes/all_ZNF_FIXED.fasta"

TSIZES="/hive/data/genomes/$GENOME/chrom.sizes"
QSIZES="/cluster/home/ifiddes/all_ZNF_FIXED.sizes"
SCORE="/scratch/data/blastz/human_chimp.v2.q"
MINSIZES="5000" #minimum size of fasta entry to save in end
MAXGAP="10000" #maximum gap size to fill in netToBed
MINSCORE="3000" #minimum chain score in axtChain

mkdir $GENOME

lastz_32 --format=axt --ambiguous=n --scores=$SCORE \
$TARGET[multiple] $QUERY[multiple] > $GENOME/$GENOME.axt

#chain axt
axtChain -verbose=0 -scoreScheme=$SCORE  \
-minScore=$MINSCORE -linearGap=medium $GENOME/$GENOME.axt $TARGET $QUERY \
stdout | chainAntiRepeat $TARGET $QUERY stdin $GENOME/$GENOME.chains

#make nets
chainPreNet $GENOME/$GENOME.chains $TSIZES $QSIZES stdout | \
chainNet stdin -minSpace=1 $TSIZES $QSIZES stdout /dev/null | \
netSyntenic stdin $GENOME/$GENOME.noClass.net

#turn net into bed, filter for fasta size
echo "track name=$GENOME description=\"$GENOME best nets\"" \
> $GENOME/$GENOME.bed
netToBed -maxGap=$MAXGAP $GENOME/$GENOME.noClass.net stdout | python \
/cluster/home/ifiddes/code/bed_region_size.py $MINSIZES >> \
$GENOME/$GENOME.bed

#make fasta from bed, rename fasta records so we know regions
fastaFromBed -s -fi $TARGET_FASTA -fo $GENOME/$GENOME.fasta -bed \
$GENOME/$GENOME.bed
python /cluster/home/ifiddes/code/fix_fasta_names_region_size.py \
$GENOME/$GENOME.fasta $GENOME/$GENOME.bed > $GENOME/$GENOME.fixed.fasta
mv $GENOME/$GENOME.fixed.fasta $GENOME/$GENOME.fasta

