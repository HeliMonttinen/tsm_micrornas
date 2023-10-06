
# This script goes through all the files in a given directory.
# and runs fpa each node pair sequences in the given files.


DIR=./results/mouse_mirbasedb2/fpa_results/
DIR3=./results/mouse_mirbasedb2/alignments2/
DIR2=./results/mouse_mirbasedb2/Ensembl_seq/*
DIR4=./results/mouse_mirbasedb2/alignments2/


declare -A pids=( )

num_procs=60


for i in $DIR2; do

	B="$(echo $i | cut -d'/' -f5)"

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	if [ ! -e "$DIR""$B" ]; then
		echo $B

		python3 ./scripts/python/main/run_fpa.py "$B" "$DIR" "$DIR3" "$DIR4" & pids["$!"]=1

	fi
done
