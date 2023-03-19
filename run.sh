#!/bin/bash

data=("email-Enron.txt" "email-Eu-core.txt" "facebook_combined.txt" "u10m_80m.txt" "musae_facebook_edges.csv")

mpic++ main.cc

for graph in ${data[@]}
do
    echo "Starting ${graph//\.txt/}"
    mpirun --use-hwthread-cpus -c 8 a.out ../data/${graph}
done
