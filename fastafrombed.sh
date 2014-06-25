for bed in *bed; do
	arr=(${bed//./ })
	species=${arr[0]}
	if [[ ! -f ~/fastas/$species.fa ]]; then
		twoBitToFa /cluster/data/$species/$species.2bit ~/fastas/$species.fa
	fi
	fastaFromBed -fi ~/fastas/$species.fa -fo $species.fasta -bed $bed -name
done
