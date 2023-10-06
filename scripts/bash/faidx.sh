
#Index genome files with samtools.

DIR=data/genomes/
DIR2=data/genomes/*fasta.gz


declare -A pids=( )

num_procs=10


for i in $DIR2; do

	echo $i


	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	samtools faidx $i & pids["$!"]=1
	
done
