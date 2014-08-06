export PASTA_DEBUG=TRUE
export PASTA_LOGGING_LEVEL=DEBUG
export PASTA_LOGGING_FORMAT=RICH

while read f; do
    name=`basename $f | cut -d "." -f 1`
    treepath="trees/$name/"
    python /cluster/home/ifiddes/pasta_1/pasta/run_pasta.py --max-mem-mb=500000 -i $f --num-cpus=18 -o $treepath -j $name --iter-limit=5 --time-limit=43200 2>$name.err >$name.log
done < $1
