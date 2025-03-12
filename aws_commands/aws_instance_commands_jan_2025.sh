
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



#building docker images
docker build -t decipherc2c-docker:1.0.5 -f Dockerfile_decipherc2c_docker_1.0.5 .
docker build -t manuscript_pre_processing:1.0.3 -f Dockerfile_manuscript_pre_processing .

#pre-processing
# Where we store the renv cache *on the host* 
# (could be any host path of your choice)
export RENV_PATHS_CACHE_HOST=/opt/local/renv/cache
# The path *inside* the container that we will mount it to.
export RENV_PATHS_CACHE_CONTAINER=/renv/cache
docker run -it --rm \
    --memory=56g --memory-swap=60g \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -v "$(pwd):/app" \
    -w /app \
    manuscript_pre_processing:1.0.3 \
    bash
    
source("scripts/")

#running on a c6a.4xlarge AWS EC2 instance
#takes about three hours to run all the stuff below, it's key to split the cpu cores efficiently as 
#if two processes or more are given the access to the same cores they may bottleneck each other due to 
#PCA calculation by CellOracle.
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
taskset -c 0-3 python3 scripts/BCG/bcg_1_cell_oracle_feb_2025.py #done
taskset -c 4-7 python3 scripts/erp/erp_1_cell_oracle_feb_2025.py #done
taskset -c 8-11 python3 scripts/tnbc/tnbc_1_cell_oracle_feb_2025.py #done
taskset -c 11-14 python3 scripts/sepsis/sepsis_1_cell_oracle_feb_2025.py #done
taskset -c 19-22 python3 scripts/covid/covid_1_cell_oracle_feb_2025.py #done
#CellOracle with more cores for PCA and KNN due to number of cells
docker run -it -m 38g --memory-swap 41g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
taskset -c 10-18 python3 scripts/lupus/lupus_1_cell_oracle_feb_2025.py #done, almost three hours

#Decipher
sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5
source("scripts/BCG/bcg_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/erp/erp_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/tnbc/tnbc_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/sepsis/sepsis_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/cord_pic/cord_pic_2_decipher_pipeline_v1_modularized_jan_2025.R") 
source("scripts/5yr_pic/5yr_pic_decipher_pipeline_v1_modularized_jan_2025.R") 
docker run -it -m 40g --memory-swap 44g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5
source("scripts/lupus/lupus_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/covid/covid_2_decipher_pipeline_v1_modularized_feb_2025.R") 

docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace ebasto/connectome:latest
source("scripts/5yr_pic/5yr_pic_3_connectome_analysis_feb_2025.R") #done
source("scripts/cord_pic/cord_pic_3_connectome_analysis_feb_2025.R") #done
source("scripts/erp/erp_3_connectome_analysis_feb_2025.R") #done but requires 60 gb of RAM
source("scripts/sepsis/sepsis_3_connectome_analysis_feb_2025.R") #done
source("scripts/BCG/bcg_3_connectome_analysis_feb_2025.R") #done
source("scripts/covid/covid_3_connectome_analysis_feb_2025.R") #done
source("scripts/lupus/lupus_3_connectome_analysis_feb_2025.R") #done but requires 130 gb or RAM
source("scripts/tnbc/tnbc_3_connectome_analysis_feb_2025.R") #done

#note that NATMI with docker requires downloading the databases and python scripts from the github page and placing them in a folder in the analysis directory called NATMI
#github https://github.com/asrhou/NATMI
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/5yr_pic/5yr_pic_5_natmi_analysis_feb_2025.sh #done
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/tnbc/tnbc_5_natmi_analysis_feb_2025.sh #done
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/sepsis/sepsis_5_natmi_analysis_feb_2025.sh #done
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/covid/covid_5_natmi_analysis_feb_2025.sh #done
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/cord_pic/cord_pic_5_natmi_analysis_feb_2025.sh #done
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/BCG/bcg_5_natmi_analysis_feb_2025.sh #done
docker run -it -m 40g --memory-swap 45g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/erp/erp_5_natmi_analysis_feb_2025.sh #done

docker run -it -m 50g --memory-swap 55g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash scripts/lupus/lupus_5_natmi_analysis_feb_2025.sh #running

