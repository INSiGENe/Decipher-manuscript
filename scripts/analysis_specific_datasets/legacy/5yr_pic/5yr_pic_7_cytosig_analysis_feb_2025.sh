#!/bin/bash

# Define the base directory
base_dir="results/5yr_pic/cytosig"

# Loop through each directory in the base directory
for dir in ${base_dir}/*/
do
    # Remove trailing slash from directory name
    dir=${dir%*/}
    
    # Get the folder name
    folder_name=$(basename ${dir})
    
    # Define input and output paths
    input_path="${dir}/differential_profile.tsv.gz"
    output_path="${dir}/output"
    
    # Run the command
    CytoSig_run.py -i ${input_path} -o ${output_path} -e 1
done

python3 scripts/silent/5yr_pic_cytosig_process_files_feb_2025.py 