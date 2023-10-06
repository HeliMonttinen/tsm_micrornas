
# This script goes through all the files in a given directory.
# and runs fpa for trees.


declare -A pids=( )

num_procs=10

for i in {1..100}; do

	while (( ${#pids[@]} >= num_procs )); do
                wait -n
                for pid in "${!pids[@]}"; do
                        kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
                done
        done

	python3 scripts/python/main/background_coordinates_bedtools.py results/genome_wide_background_corrected_tsms/background $i results/unmasked_variant_reg_tsms.txt all.mask2.bed & pids["$!"]=1

done
