"""
Parses identifiers from the Emseml_seq
directory and aligns the corresponding
files in the results/mouse_mirbasedb2/alignments2/ directory
through the script ./scripts/python/main/alignments2.py
"""
DIR=./results/mouse_mirbasedb2/alignments2/
DIR2=./results/mouse_mirbasedb2/Ensembl_seq/*



declare -A pids=( )

num_procs=20


for i in $DIR2; do

	B="$(echo $i | cut -d'/' -f5)"
	#echo $B"_pagan.fas"

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	if [ ! -e "$DIR""$B""_pagan.fas" ]; then
		python3 ./scripts/python/main/alignments2.py "$DIR" "$DIR" "$DIR" "$B" & pids["$!"]=1
	
	fi
done
