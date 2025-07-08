
#################################
####### Download Docker images ############
#################################
#TODO: add docker images to dockerHub and include commands to download docker files here

####### Alternative: Build Docker images ############
docker build -t decipherc2c-docker:1.0.5 -f Dockerfile_decipherc2c_docker_1.0.5 .
docker build -t manuscript_pre_processing:1.0.3 -f Dockerfile_manuscript_pre_processing .

#################################
####### Analysis ############
#################################
#navigate to analysis directory
mkdir projects/analysis
cd projects/analysis

# annotate datasets without annotations
# run azimuth on SevMildCOVID
cd data/SevMildCOVID
docker run -it --rm -v "$(pwd):/workspace" -w /workspace satijalab/azimuth:0.5.0 bash

#### Convert anndata objects to R-based objects (Seurat) ####
#trigger environment
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest

#general command form
python3 scripts/1_preprocess_h5ad.py dataset_key

#commands for each dataset
python3 scripts/1_preprocess_h5ad.py lupus

#### Process initial object for downstream analyses ####
#TODO: convert this to a static docker image (now that analysis is complete)
export RENV_PATHS_CACHE_HOST=/opt/local/renv/cache
# The path *inside* the container that we will mount it to.
export RENV_PATHS_CACHE_CONTAINER=/renv/cache
docker run -it --rm \
    --memory=180g --memory-swap=185g \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -v "$(pwd):/app" \
    -w /app \
    manuscript_pre_processing:1.0.3 \
    bash

#### --------------------------------- ####
####  Custom pre-processing (selected) ####
#### --------------------------------- ####
# run azimuth on SevMildCOVID
#TODO: understand why we have three different SevMilCOVID folders here?
Rscript scripts/custom_pre_processing_SevMilCovid.r SevCOVID
Rscript scripts/custom_pre_processing_SevMilCovid.r MilCOVID

Rscript scripts/1.1_preprocess_SevMilCovid_Azimuth.r SevCOVID_Azimuthl1
Rscript scripts/1.1_preprocess_SevMilCovid_Azimuth.r MilCOVID_Azimuthl1

Rscript scripts/1.1_preprocess_SevMilCovid_Azimuth.r SevCOVID_Azimuthl2
Rscript scripts/1.1_preprocess_SevMilCovid_Azimuth.r MilCOVID_Azimuthl2

#### ----------------------------- ####
####  Generic pre-processing (all) ####
#### ----------------------------- ####
#TODO: distinguish between CZ pipeline and custom pipelines

Rscript scripts/2_preprocess_object_for_analysis.R dataset_key

#### run scCODA analysis for SevMildCOVID ####
docker run -it \
  -v "$(pwd):/workspace" \
  -w /workspace \
  wollmilchsau/scanpy_sccoda:latest
python3
#TODO: check this python3 command here

#TODO: operate from a single-directory, rather than move folders
#TODO: here what's being moved is the scoda output I guess? check
sudo mkdir -p Manuscript_jan_2025/results/SevCOVID_Azimuthl2/sccoda
sudo mkdir -p Manuscript_jan_2025/results/MilCOVID_Azimuthl2/sccoda
sudo mv pre_processing_test/data/SevMilCOVID/results_sccoda_severe_vs_healthy.csv Manuscript_jan_2025/results/SevCOVID_Azimuthl2/sccoda
sudo mv pre_processing_test/data/SevMilCOVID/results_sccoda_moderate_vs_healthy.csv Manuscript_jan_2025/results/MilCOVID_Azimuthl2/sccoda

#### Move results to analysis folder ####
#TODO: operate from a single-directory, rather than move folders
#generic command
sudo mv pre_processing_test/results/dataset_key Manuscript_jan_2025/results/

#commands for all datasets
sudo mv pre_processing_test/results/MilCOVID_Azimuthl2 Manuscript_jan_2025/results/


#### ----------------------- ####
####       Run Cytosig       ####
#### ----------------------- ####
#TODO: operate from a single analysis directory
cd projects/Manuscript_jan_2025

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite

#generic analysis command
bash scripts/3_cytosig_run.sh dataset_key

#commands for all datasets
bash scripts/3_cytosig_run.sh MilCOVID_Azimuthl2

#### ----------------------- ####
#### Run CellOracle analysis ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest

#generic analysis command
#NOTE: if you want to run multiple CellOracle analyses in parallel, constraining the number of cores per analysis prevents overflow
taskset -c 0-3 python3 scripts/5_cell_oracle.py dataset_key

#example of running two datasets in parallel (initialize a distinct docker container for each)
taskset -c 0-3 python3 scripts/5_cell_oracle.py SevCOVID_Azimuthl2
taskset -c 4-7 python3 scripts/5_cell_oracle.py MilCOVID_Azimuthl2

#commands for all datasets

#### ----------------------- ####
####      Run Decipher       ####
#### ----------------------- ####
#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5 bash

#generic analysis command
Rscript scripts/6_decipher_pipeline_v1_modularized.R dataset_key

#commands for all datasets
Rscript scripts/6_decipher_pipeline_v1_modularized.R SevCOVID_Azimuthl2_k0


#### ----------------------- ####
####      Run Connectome     ####
#### ----------------------- ####
#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/connectome:latest bash

#generic analysis command
Rscript scripts/7_connectome_analysis.R dataset_key

#commands for all datasets

#### ----------------------- ####
####      Run NicheNet       ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest bash

#generic analysis command
Rscript scripts/8_nichenet_analysis.R dataset_key

#### ----------------------- ####
####        Run NATMI        ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace asrhou/natmi  

#generic analysis command
bash scripts/9_natmi_analysis.sh dataset_key

#### ----------------------- ####
####        Run LIANA+       ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus bash

#generic analysis command
python3 scripts/10_liana_plus_analysis.py dataset_key