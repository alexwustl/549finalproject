for bytes in 64 1024
do
for i in 9 25 49
do
    res=""
    res2=""
    for j in 100000 500000 1000000 5000000 10000000
    do
        if [[ $i -eq 49 ]]
        then
            k=7
        elif [[ $i -eq 25 ]]
        then
                k=5
        else
                k=3
        fi
        gcc -DVALUE_BYTES=$bytes -DNUM_ITERS=$j -DNUM=$i -DBEPSILON=$k red_black3.c btreeepsilon.c
        perf stat --repeat=3 -B -e cache-references,cache-misses,L1-dcache-load-misses,LLC -o temp.txt ./a.out
        restemp=$(cat temp.txt | grep -oP "[0-9]+(?= \s+cache-misses)")
        res="${res}, ${restemp}"
        gcc -DVALUE_BYTES=$bytes -DNUM_ITERS=$j -DNUM=$i btree.c
        perf stat --repeat=3 -B -e cache-references,cache-misses,L1-dcache-load-misses,LLC -o temp.txt ./a.out
        restemp=$(cat temp.txt | grep -oP "[0-9]+(?= \s+cache-misses)")
        res2="${res2}, ${restemp}"
    done
    echo "Epsilon" $res
    echo "BTree" $res2
done
done
