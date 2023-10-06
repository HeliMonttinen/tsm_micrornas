
# This script goes through all the files in a given directory.
# and runs fpa each node pair sequences in the given files.


DIR=./results/fpa_microRNA/
DIR3=./results/alignments2/
DIR2=./results/mafseqs_full/*
DIR4=./results/alignments2/


declare -A pids=( )

num_procs=60


for i in $DIR2; do

	B="$(echo $i | cut -d'/' -f4)"

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
