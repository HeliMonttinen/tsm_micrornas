
# gzip to bgzip.


DIR=data/genomes/
DIR2=data/genomes/*



declare -A pids=( )

num_procs=10


for i in $DIR2; do

	echo $i

	A="$(echo $i | cut -d'.' -f1)"
	echo $A

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	if [ ! -e "$A"".fasta.gz" ]; then

		zcat  "$i" | bgzip -c > "$A"".fasta.gz" & pids["$!"]=1

	fi
done
