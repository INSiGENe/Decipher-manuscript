#set the working directory to the project root (Decipher folder)

#################################
####### Download Docker images ############
#################################
docker pull ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a
docker pull ebasto/decipher-manuscript-cytosig@sha256:583a450ae25f91686dbf9db9b3297d7f82f81f87b809f520daf226c3a661b11b
docker pull ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1
docker pull ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5
docker pull ebasto/decipher-manuscript-connectome@sha256:851edafd24ec10d67cd9eda8c07f611a22baf13da5667acd539bb8a390c05973
docker pull ebasto/decipher-manuscript-nichenetr@sha256:146e50752019a18d4125729c3d4b09ccd88a11e188525dc48eb9dcdcb72967ec
docker pull ebasto/decipher-manuscript-natmi@sha256:1269438fa798330eba47d51ac910d76a6298fb471c6e9449685a0a1dbb2282b7 
docker pull ebasto/decipher-manuscript-liana-plus:1.0.0@sha256:d300ec7872d9a0cf8ae91fc05798f56a3aa3982657bb3883f21bcef63b8ee580

#################################
#######  Manuscript pre-processing ############
#################################
sudo docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a
source("scripts/analysis_specific_datasets/5yr_pic/5yr_pic_0_pre_processing_feb_2025.R") 
source("scripts/analysis_specific_datasets/BCG/bcg_0_pre_processing_feb_2025.R") 
source("scripts/analysis_specific_datasets/cord_pic/cord_pic_0_pre_processing_jan_2025.R") 
source("scripts/analysis_specific_datasets/covid/covid_0_pre_processing_feb_2025.R") 
source("scripts/analysis_specific_datasets/erp/erp_0_pre_processing_feb_2025.R") 
source("scripts/analysis_specific_datasets/lupus/lupus_0_pre_processing_feb_2025.R") 
source("scripts/analysis_specific_datasets/sepsis/sepsis_0_pre_processing_feb_2025.R") 
source("scripts/analysis_specific_datasets/tnbc/tnbc_0_pre_processing_feb_2025.R") 

#################################
####### CellOracle ############
#################################
#running on a c6a.4xlarge AWS EC2 instance
#takes about three hours to run all the stuff below, it's key to split the cpu cores efficiently as 
#if two processes or more are given the access to the same cores they may bottleneck each other due to 
#PCA calculation by CellOracle.
docker run -it -m 60g --memory-swap 62g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1
taskset -c 0-3 python3 scripts/analysis_specific_datasets/BCG/bcg_1_cell_oracle_feb_2025.py
taskset -c 4-7 python3 scripts/analysis_specific_datasets/erp/erp_1_cell_oracle_feb_2025.py
taskset -c 8-11 python3 scripts/analysis_specific_datasets/tnbc/tnbc_1_cell_oracle_feb_2025.py 
taskset -c 11-14 python3 scripts/analysis_specific_datasets/sepsis/sepsis_1_cell_oracle_feb_2025.py 
taskset -c 19-22 python3 scripts/analysis_specific_datasets/covid/covid_1_cell_oracle_feb_2025.py 
tasket -c 23-25 python3 scripts/analysis_specific_datasets/5yr_pic/5yr_pic_1_cell_oracle_jan_2025.py
tasket -c 26-28 python3 scripts/analysis_specific_datasets/cord_pic/cord_pic_1_cell_oracle_jan_2025.py
taskset -c 10-18 python3 scripts/analysis_specific_datasets/lupus/lupus_1_cell_oracle_feb_2025.py #this requires more cores due to the size of the dataset

#################################
####### Decipher ############
#################################
sudo docker run -it -m 40g --memory-swap 44g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5
#select yes when asked about setting up reticulate and no if asked about muscData

source("scripts/analysis_specific_datasets/BCG/bcg_2_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/erp/erp_2_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/tnbc/tnbc_2_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/sepsis/sepsis_2_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/cord_pic/cord_pic_2_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/5yr_pic/5yr_pic_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/lupus/lupus_2_decipher_pipeline_v1_modularized.R") 
source("scripts/analysis_specific_datasets/covid/covid_2_decipher_pipeline_v1_modularized.R") 

#################################
####### Connectome ############
#################################
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-connectome@sha256:851edafd24ec10d67cd9eda8c07f611a22baf13da5667acd539bb8a390c05973
source("scripts/analysis_specific_datasets/5yr_pic/5yr_pic_3_connectome_analysis_feb_2025.R") 
source("scripts/analysis_specific_datasets/cord_pic/cord_pic_3_connectome_analysis_feb_2025.R") 
source("scripts/analysis_specific_datasets/erp/erp_3_connectome_analysis_feb_2025.R") #requires 60 gb of RAM
source("scripts/analysis_specific_datasets/sepsis/sepsis_3_connectome_analysis_feb_2025.R") 
source("scripts/analysis_specific_datasets/BCG/bcg_3_connectome_analysis_feb_2025.R") 
source("scripts/analysis_specific_datasets/covid/covid_3_connectome_analysis_feb_2025.R") 
source("scripts/analysis_specific_datasets/lupus/lupus_3_connectome_analysis_feb_2025.R") #requires 130 gb or RAM
source("scripts/analysis_specific_datasets/tnbc/tnbc_3_connectome_analysis_feb_2025.R") 

#################################
####### NATMI ############
#################################
#note that NATMI with docker requires downloading the databases and python scripts from the github page and placing them in a folder in the analysis directory called NATMI
#github https://github.com/asrhou/NATMI
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-natmi@sha256:1269438fa798330eba47d51ac910d76a6298fb471c6e9449685a0a1dbb2282b7 
bash scripts/analysis_specific_datasets/5yr_pic/5yr_pic_5_natmi_analysis_feb_2025.sh
bash scripts/analysis_specific_datasets/tnbc/tnbc_5_natmi_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/sepsis/sepsis_5_natmi_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/covid/covid_5_natmi_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/cord_pic/cord_pic_5_natmi_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/BCG/bcg_5_natmi_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/erp/erp_5_natmi_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/lupus/lupus_5_natmi_analysis_feb_2025.sh #running

#################################
####### NicheNet ############
#################################
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-nichenetr@sha256:146e50752019a18d4125729c3d4b09ccd88a11e188525dc48eb9dcdcb72967ec
Rscript scripts/analysis_specific_datasets/5yr_pic/5yr_pic_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/BCG/bcg_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/cord_pic/cord_pic_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/covid/covid_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/erp/erp_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/tnbc/tnbc_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/sepsis/sepsis_4_nichenet_analysis_feb_2025.R 
Rscript scripts/analysis_specific_datasets/lupus/lupus_4_nichenet_analysis_feb_2025.R 

#################################
####### LIANA+ ############
#################################

docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-liana-plus:1.0.0@sha256:d300ec7872d9a0cf8ae91fc05798f56a3aa3982657bb3883f21bcef63b8ee580 
python3 scripts/analysis_specific_datasets/5yr_pic/5yr_pic_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/BCG/bcg_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/cord_pic/cord_pic_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/covid/covid_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/erp/erp_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/tnbc/tnbc_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/sepsis/sepsis_6_liana_plus_analysis_feb_2025.py 
python3 scripts/analysis_specific_datasets/lupus/lupus_6_liana_plus_analysis_feb_2025.py 


#################################
####### Cytosig run and process data ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-cytosig@sha256:583a450ae25f91686dbf9db9b3297d7f82f81f87b809f520daf226c3a661b11b
bash scripts/analysis_specific_datasets/5yr_pic/5yr_pic_7_cytosig_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/BCG/bcg_7_cytosig_analysis_feb_2025.sh #running
bash scripts/analysis_specific_datasets/cord_pic/cord_pic_7_cytosig_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/covid/covid_7_cytosig_analysis_feb_2025.sh #running
bash scripts/analysis_specific_datasets/erp/erp_7_cytosig_analysis_feb_2025.sh 
bash scripts/analysis_specific_datasets/lupus/lupus_7_cytosig_analysis_feb_2025.sh #running
bash scripts/analysis_specific_datasets/sepsis/sepsis_7_cytosig_analysis_feb_2025.sh #running
bash scripts/analysis_specific_datasets/tnbc/tnbc_7_cytosig_analysis_feb_2025.sh



#################################
####### Robustness benchmarking ############
#################################
#download and create sample seurat object
#select yes when asked about setting up reticulate and no if asked about muscData
sudo docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5
source('scripts/analysis_specific_datasets/sample_dataset/sample_dataset_GRN.R') 

#prepare data for CellOracle
docker run -it --rm -v "$(pwd):/app" -w /app ebasto/manuscript_pre_processing:1.0.4@sha256:9b5c93bba509359a11181bbb297e1af3b99c8b3130e8adb105549414ebd0fb0a bash
Rscript scripts/analysis_specific_datasets/sample_dataset/write_data_for_CellOracle.r

#run CellOracle
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1
taskset -c 0-3 python3 scripts/analysis_cellxgene_datasets/5_cell_oracle.py sample_analysis

#run benchmarking
sudo docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5
source('scripts/benchmarking_and_figure_scripts/2a_robustness_benchmarking.R') 

