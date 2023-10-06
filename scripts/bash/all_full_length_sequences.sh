#!/bin/bash
# This script goes through all the files in a given directory.
# and retrieves full lengths sequences for each chunk.


DIR=./results/maf_sequence_all/*fas


declare -A pids=( )

num_procs=30

for i in $DIR; do
	

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done


	
	python3 scripts/python/main/ucsc_seq_all.py data/species_mammals.txt results/sequence_info_all $i results/maf_sequence_full_all & pids["$!"]=1

done
