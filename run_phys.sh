#! /bin/bash

rm -rf results/*

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ariasn/gsl/lib && export LD_LIBRARY_PATH
NUM_RUNS=10.0
SUM=0
for i in $(seq 1 10)
do
    RUN=$(./sequential)
    SUM=$(echo "$RUN + $SUM" | bc)
done
echo "Sequential, $(echo "scale=4; $SUM/$NUM_RUNS" | bc)"

for PROCESS in $(seq 1 112)
do
    SUM=0
    for i in $(seq 1 10)
    do
        RUN=$(./run_parallel.sh $PROCESS)
        SUM=$(echo "$RUN + $SUM" | bc)
    done
    AVG=$(echo "scale=4; $SUM/$NUM_RUNS" | bc)
    echo "$PROCESS, $AVG"
done
#paste -d , results/* > consolidated_results.csv
