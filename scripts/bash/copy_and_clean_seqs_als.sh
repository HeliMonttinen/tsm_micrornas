
# This script goes through all the files in a given directory.
# the script identifies very different sequences and copies them
# to another directory.


DIR=./pubset/results/alignments1_all/
DIR4=./pubset/results/alignments2_all/
DIR2=./pubset/results/maf_sequence_full_all/*
DIR3=./pubset/results/maf_sequence_full_all/


declare -A pids=( )

num_procs=30


for i in $DIR2; do

	A="$(echo $i | cut -d'.' -f1)"
	B="$(echo $A | cut -d'/' -f8)"
	#echo $B"_pagan.fas"

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	python3 ./scripts/python/main/copy_seqs_all_genes.py "$DIR4" "$DIR" "$B" & pids["$!"]=1


done
