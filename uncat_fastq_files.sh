#!/usr/bin/bash

for f in *; do
    len=`wc -l $f | cut -d " " -f 1`
    len=$((len/2))
    name=`echo $f | cut -d "." -f 1`
    head -n $len $f > $name.R1.fastq
    tail -n $len $f > $name.R2.fastq
done