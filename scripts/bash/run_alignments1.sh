
# This script goes through all the files in a given directory.
# and runs pagan alignments for each file if _pagan.fas file
# does not exist yet.


DIR=./results/alignments1/
DIR2=./results/mafseqs_full/*
DIR3=./results/alignments1/



declare -A pids=( )

num_procs=60


for i in $DIR2; do

	B="$(echo $i | cut -d'/' -f4)"
	echo $B

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done
	echo $DIR$B
	if [ ! -e "$DIR""$B""_pagan.fas" ]; then
		python3 ./scripts/python/main/alignments1.py "$DIR" "$DIR" "$DIR" "$B" & pids["$!"]=1

	fi
done
