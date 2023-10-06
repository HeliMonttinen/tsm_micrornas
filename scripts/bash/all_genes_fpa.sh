# This script goes through all the files in a given directory
# and runs fpa for each file.


DIR=./results/fpa_all/
DIR3=./results/alignments2_all/
DIR2=./results/maf_sequence_full_all/*
DIR4=./results/alignments2_all/


declare -A pids=( )

num_procs=30


for i in $DIR2; do

	B="$(echo $i | cut -d'/' -f4)"
	echo $B

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	if [ ! -e "$DIR""$B" ]; then

		python3 ./scripts/python/main/run_fpa.py "$B" "$DIR" "$DIR3" "$DIR4" & pids["$!"]=1

	fi
done
