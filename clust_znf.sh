#!/bin/bash
######################################################
#Run lastz on parasol on two directories of single-entry fastas
#Expects there to be two directories in the folder with single-entry fastas
#representing each of the zinc fingers and each entry in a genome
#use split_fa_to_seqs to split up genome fastas
######################################################


#check for 2 input arguments
if [ $# -ne 2 ]; then
	echo "Two arguments required. 1st - directory of fastas (big one - genome). 2nd - directory of fastas (small one - ZNF)"
	exit 1
fi

#variables to change job directory/combined job easier
JOBDIR=jobDir
CATJOB=jobList

#build array of genome fasta file names
GENOME=( ` find $1/ -name '*.fa' -o -name '*.fasta' ` )

#make sure $JOBDIR does not exist; otherwise may have overlaps
if [ -d $JOBDIR ]; then
	echo "Error: $JOBDIR exists. Exiting"
	exit 1
fi

mkdir $JOBDIR


#iterate over ZNF fasta folder
for FASTA in $2/*; do
		#do sanity checks; since we are using >> to concatenate
        FASTBASE=`basename $FASTA`
        if [ -f $FASTBASE.jobList ]; then
        	echo "Error: $FASTBASE.jobList exists. Exiting"
        	exit 1
        fi
        if [ -d $FASTBASE ]; then
        	echo "Error: $FASTBASE exists. Exiting"
        	exit 1
        fi
        mkdir $FASTBASE
        touch $JOBDIR/$FASTBASE.jobList
        for CHR in ${GENOME[@]}; do
        		CHRBASE=`basename $CHR`
                echo "lastz_32 --format=general:name1,zstart1,end1,name2,zstart2+,end2+,size2,identity,coverage,length1,length2,nmatch,nmismatch,cigarx- --chain --ambiguous=n --output=$FASTBASE/$FASTBASE.$CHRBASE.result $FASTA $CHR" >> $JOBDIR/$FASTBASE.jobList
        done
done

#run on parasol
echo '#!/bin/csh -fe' > $CATJOB
cat $JOBDIR/*.jobList >> $CATJOB
para create $CATJOB
para shove

#combine all results together, parse for empty files
#python2.7 /cluster/home/ifiddes/code/combine_clust.py > combined.results


