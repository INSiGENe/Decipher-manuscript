sudo docker build -t decipherc2c-docker:1.0.3 -f Dockerfile_w_python_v2 .


#activate celloracle env in a net terminal window and run
#change to analysis directory
python3 scripts/sepsis_1_cell_oracle.py

#activate natmi env
bash scripts/sepsis_5_natmi_analysis.sh

#activate liana env
python3 scripts/sepsis_6_liana_plus_analysis.py



#for cytosig
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite
bash scripts/sepsis_7_cytosig_analysis.sh
python3 scripts/silent/sepsis_cytosig_process_files.py


docker run -it -v "$(pwd):/workspace" -w /workspace kenjikamimoto126/celloracle_ubuntu:0.18.0

docker run -it -m 12g --memory-swap 16g -v "$(pwd):/workspace" -w /workspace kenjikamimoto126/celloracle_ubuntu:0.18.0




#for celloracle

docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest




#5YR PIC

sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.3
source("scripts/5yr_pic/5yr_pic_0_pre_processing_jan_2025.R")

docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
python3 /scripts/5yr_pic/5yr_pic_1_cell_oracle_jan_2025.py

sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.3
source("scripts/5yr_pic/5yr_pic_2_decipher_pipeline_v1_modularized_jan_2025.R")

#CORD PIC

sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.3
source("scripts/cord_pic/cord_pic_0_pre_processing_jan_2025.R")

docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
python3 /scripts/cord_pic/cord_pic_1_cell_oracle_jan_2025.py

sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.3
source("scripts/cord_pic/cord_pic_2_decipher_pipeline_v1_modularized_jan_2025.R")

#BCG

sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.4
source("scripts/BCG/bcg_0_pre_processing_jan_2025.R")


docker build -t manuscript_pre_processing:1.0.2 -f Dockerfile_manuscript_pre_processing .
sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace r-base:4.4.1

sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace --entrypoint /bin/bash manuscript_pre_processing:1.0.1






# Where we store the renv cache *on the host* 
# (could be any host path of your choice)
export RENV_PATHS_CACHE_HOST=/opt/local/renv/cache

# The path *inside* the container that we will mount it to.
export RENV_PATHS_CACHE_CONTAINER=/renv/cache

docker run -it --rm \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -v "$(pwd):/app" \
    -w /app \
    manuscript_pre_processing:1.0.2
    
    
     \
    -e 'renv::restore(); print("All packages restored!"); # your code here'
