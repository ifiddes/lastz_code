#!/bin/bash

read -ra r1 <<< `find $1 | grep "R1"`
read -ra r2 <<< `find $1 | grep "R2"`

for ((i=0;i<${#r1[@]};i++));do
	name=`echo ${r1[i]} | grep -o '[^/]*$' | cut -d "_" -f 1 | cut -d "-" -f 2`
	pear -f ${r1[i]} -r ${r2[i]} -o $name
done
