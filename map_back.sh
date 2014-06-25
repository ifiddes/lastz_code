for bed in *bed; do
	arr=(${bed//./ })
	species=${arr[0]}
	if [[ ! -f ~/fastas/$species.fa ]]; then
		twoBitToFa /cluster/data/$species/$species.2bit ~/fastas/$species.fa
	fi
	fastaFromBed -fi ~/fastas/$species.fa -fo $species.fasta -bed $bed -name
done

for fasta in *fasta; do
	arr=(${fasta//./ })
	species=${arr[0]}
	python2.7 ~/code/lastz_32_general_format.py --target /cluster/data/hg38/hg38.2bit --cores 16 --chain --fasta $fasta > $species.hg38.mapback.chains
	python2.7 ~/code/merge_chains.py $species.hg38.mapback.chains > $species.hg38.mapback.chains.merged
	python2.7 ~/code/merged_to_bed.py $species.hg38.mapback.chains.merged > $$species.hg38.mapback.chains.merged.bed
done
