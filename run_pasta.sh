while read f; do
    name=`basename $f | cut -d "." -f 1`
    treepath="trees/$name/"
    python /cluster/home/ifiddes/pasta_1/pasta/run_pasta.py --max-mem-mb=500000 -i $f --num-cpus=18 -o $treepath -j $name --time-limit=43200
done < $1
