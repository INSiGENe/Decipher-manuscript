#!/bin/bash

# Usage: ./run_cytosig.sh cz_rcc
set -e

# Dataset key passed as argument
dataset_key=$1

# Extract output path using jq
output_path=$(python3 -c "
import json, sys
with open('scripts/config.json') as f:
    data = json.load(f)
    print(data.get('$dataset_key', {}).get('cytosig', {}).get('output_path', ''))
")

if [ -z "$output_path" ] || [ "$output_path" == "null" ]; then
    echo "❌ No output_path found for dataset key: ${dataset_key}"
    exit 1
fi

# Construct base directory for cytosig
base_dir="${output_path}/cytosig"

# Loop through each cluster folder
for dir in ${base_dir}/*/
do
    dir=${dir%*/}  # Remove trailing slash
    folder_name=$(basename "$dir")

    input_path="${dir}/differential_profile.tsv.gz"
    output_path="${dir}/output"

    echo "Running CytoSig on $folder_name"
    CytoSig_run.py -i "$input_path" -o "$output_path" -e 1
done

# Run follow-up processing step
python3 scripts/cytosig_clean.py dataset_key
