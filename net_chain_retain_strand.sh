#!/bin/bash


GENOME=`basename $1 | cut -d "." -f 1`
TARGET="/hive/data/genomes/$GENOME/$GENOME.2bit"
QUERY="/cluster/home/ifiddes/znf/stranded_znf.2bit"
TARGET_FASTA="/cluster/home/ifiddes/fastas/$GENOME.fa"
QUERY_FASTA="/cluster/home/ifiddes/znf/stranded_znf.fasta"

TSIZES="/hive/data/genomes/$GENOME/chrom.sizes"
QSIZES="/cluster/home/ifiddes/znf/stranded_znf.sizes"
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
netSyntenic stdin $GENOME/$GENOME.net

#hacky way to get strand information from nets - thanks Joel
#filter for $MINSIZES as well as $MAXGAP
echo "track name=$GENOME description=\"$GENOME minsize=$MINSIZES maxgap=$MAXGAP minscore=$MINSCORE\"" \
> $GENOME/$GENOME.bed

netToAxt $GENOME/$GENOME.net $GENOME/$GENOME.chains $TARGET \
$QUERY stdout | axtToBed -extended stdin stdout | \
cut -f 1,2,3,4,5,6 - | bedtools sort > /dev/stdout | \
bedtools merge -d $MAXGAP -s -nms -i /dev/stdin | \
python /cluster/home/ifiddes/code/collapse_bed_names_from_merging.py | \
python /cluster/home/ifiddes/code/bed_region_size.py $MINSIZES >> \
$GENOME/$GENOME.bed

#make fasta from bed
fastaFromBed -s -fi $TARGET_FASTA -fo $GENOME/$GENOME.fasta -bed \
$GENOME/$GENOME.bed
#rename fasta records so we know regions
python /cluster/home/ifiddes/code/fix_fasta_names_region_size.py \
$GENOME/$GENOME.fasta $GENOME/$GENOME.bed > $GENOME/$GENOME.fixed.fasta
mv $GENOME/$GENOME.fixed.fasta $GENOME/$GENOME.fasta
