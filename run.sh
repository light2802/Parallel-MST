#!/bin/bash

data=("clean-soc-pokec-relationships.txt"  "GermanyRoadud.txt" "clean-soc-sinaweibo.txt" "USAud.txt" "clean-soc-twitter.txt")

mpic++ main.cc

for graph in ${data[@]}
do
    echo "Starting ${graph//\.txt/}"
    mpirun --use-hwthread-cpus -c 8 a.out ../data/${graph}
done
