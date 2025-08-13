#set the working directory to the project root (Decipher folder)

#################################
####### Download Docker images ############
#################################
docker pull ebasto/decipher-manuscript-cytosig@sha256:583a450ae25f91686dbf9db9b3297d7f82f81f87b809f520daf226c3a661b11b
docker pull ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1
docker pull ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5
docker pull ebasto/decipher-manuscript-connectome@sha256:851edafd24ec10d67cd9eda8c07f611a22baf13da5667acd539bb8a390c05973
docker pull ebasto/decipher-manuscript-nichenetr@sha256:146e50752019a18d4125729c3d4b09ccd88a11e188525dc48eb9dcdcb72967ec
docker pull ebasto/decipher-manuscript-natmi@sha256:1269438fa798330eba47d51ac910d76a6298fb471c6e9449685a0a1dbb2282b7 
docker pull ebasto/decipher-manuscript-liana-plus:1.0.0@sha256:d300ec7872d9a0cf8ae91fc05798f56a3aa3982657bb3883f21bcef63b8ee580
docker pull ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a
docker pull ebasto/decipher-manuscript-azimuth:1.0.0@sha256:128653884fd49cf2f59b00cd2c8ca6509e6cc111784ced78f6c24b163600f4ab

#################################
####### Analysis ############
#################################

### ---> #IMPORTANT: replace dataset_key with the dataset root  of each structure  in the scripts/config.json,, e.g. cz_human_kidney_v1.5, cz_influenza, etc. 

#### --------------------------------- ####
####  pre-process python objects (All CELLxGENE except Severe vs Mild COVID_19 dataset) ####
#### --------------------------------- ####

#### Convert CELLxGENE anndata objects to Seurat objects ####
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1
#general command form
python3 scripts/analysis_cellxgene_datasets/1_preprocess_h5ad.py dataset_key

#### --------------------------------- ####
####  Custom pre-processing (Only Severe vs Mild COVID-19) ####
#### --------------------------------- ####

# run intial pre-processing of raw data
docker run -it --rm -v "$(pwd):/app" -w /app ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a bash
Rscript scripts/analysis_cellxgene_datasets/custom_sevmil_covid/1_create_seurat_object.r

# run azimuth on SevMildCOVID
docker run -it --rm -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-azimuth:1.0.0@sha256:128653884fd49cf2f59b00cd2c8ca6509e6cc111784ced78f6c24b163600f4ab bash
Rscript scripts/analysis_cellxgene_datasets/custom_sevmil_covid/2_annotation.r

# separate azimuth results into severe and mild
docker run -it --rm -v "$(pwd):/app" -w /app ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a bash
Rscript scripts/analysis_cellxgene_datasets/custom_sevmil_covid/3_separate_severe_and_mild.r SevCOVID_Azimuthl2
Rscript scripts/analysis_cellxgene_datasets/custom_sevmil_covid/3_separate_severe_and_mild.r MilCOVID_Azimuthl2

#### ----------------------------- ####
####  Generic pre-processing (all) ####
#### ----------------------------- ####
#TODO: distinguish between CZ pipeline and custom pipelines
docker run -it --rm -v "$(pwd):/app" -w /app ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a bash
  
Rscript scripts/analysis_cellxgene_datasets/2_preprocess_object_for_analysis.R dataset_key

#### ----------------------- ####
####       Run Cytosig (all)      ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-cytosig@sha256:583a450ae25f91686dbf9db9b3297d7f82f81f87b809f520daf226c3a661b11b

#generic analysis command
bash scripts/analysis_cellxgene_datasets/3_cytosig_run.sh dataset_key

#commands for all datasets
bash scripts/analysis_cellxgene_datasets/3_cytosig_run.sh MilCOVID_Azimuthl2

#### ----------------------- ####
#### Run CellOracle analysis ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1

#generic analysis command
#NOTE: if you want to run multiple CellOracle analyses in parallel, constraining the number of cores per analysis prevents overflow
taskset -c 0-3 python3 scripts/analysis_cellxgene_datasets/5_cell_oracle.py dataset_key

#example of running two datasets in parallel (initialize a distinct docker container for each)
taskset -c 0-3 python3 scripts/analysis_cellxgene_datasets/5_cell_oracle.py SevCOVID_Azimuthl2
taskset -c 4-7 python3 scripts/analysis_cellxgene_datasets/5_cell_oracle.py MilCOVID_Azimuthl2

#commands for all datasets

#### ----------------------- ####
####      Run Decipher       ####
#### ----------------------- ####
#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5 bash
#select yes when asked about setting up reticulate and no if asked about muscData

#generic analysis command
Rscript scripts/analysis_cellxgene_datasets/6_decipher_pipeline_v1_modularized.R MilCOVID_Azimuthl2


#### ----------------------- ####
####      Run Connectome     ####
#### ----------------------- ####
#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-connectome@sha256:851edafd24ec10d67cd9eda8c07f611a22baf13da5667acd539bb8a390c05973 bash

#generic analysis command
Rscript scripts/analysis_cellxgene_datasets/7_connectome_analysis.R dataset_key

#commands for all datasets

#### ----------------------- ####
####      Run NicheNet       ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-nichenetr@sha256:146e50752019a18d4125729c3d4b09ccd88a11e188525dc48eb9dcdcb72967ec bash

#generic analysis command
Rscript scripts/analysis_cellxgene_datasets/8_nichenet_analysis.R dataset_key

#### ----------------------- ####
####        Run NATMI        ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-natmi@sha256:1269438fa798330eba47d51ac910d76a6298fb471c6e9449685a0a1dbb2282b7 

#generic analysis command
bash scripts/analysis_cellxgene_datasets/9_natmi_analysis.sh dataset_key

#### ----------------------- ####
####        Run LIANA+       ####
#### ----------------------- ####

#analysis container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-liana-plus:1.0.0@sha256:d300ec7872d9a0cf8ae91fc05798f56a3aa3982657bb3883f21bcef63b8ee580 bash

#generic analysis command
python3 scripts/analysis_cellxgene_datasets/10_liana_plus_analysis.py dataset_key
