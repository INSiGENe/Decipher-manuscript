
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
#TODO: missing manuscript_pre_processing

#################################
####### Run Decipher analysis ############
#################################
sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5

#################################
####### CellOracle ############
#################################
#running on a c6a.4xlarge AWS EC2 instance
#takes about three hours to run all the stuff below, it's key to split the cpu cores efficiently as 
#if two processes or more are given the access to the same cores they may bottleneck each other due to 
#PCA calculation by CellOracle.
docker run -it -m 60g --memory-swap 62g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1
taskset -c 0-3 python3 scripts/BCG/bcg_1_cell_oracle_feb_2025.py
taskset -c 4-7 python3 scripts/erp/erp_1_cell_oracle_feb_2025.py
taskset -c 8-11 python3 scripts/tnbc/tnbc_1_cell_oracle_feb_2025.py 
taskset -c 11-14 python3 scripts/sepsis/sepsis_1_cell_oracle_feb_2025.py 
taskset -c 19-22 python3 scripts/covid/covid_1_cell_oracle_feb_2025.py 
tasket -c 23-25 python3 scripts/5yr_pic/5yr_pic_1_cell_oracle_jan_2025.py
tasket -c 26-28 python3 scripts/analysis_specific_datasets/cord_pic/cord_pic_1_cell_oracle_jan_2025.py
taskset -c 10-18 python3 scripts/lupus/lupus_1_cell_oracle_feb_2025.py #this requires more cores due to the size of the dataset

#taskset -c 0-3 python cell_oracle_mar_2025.py SkinAtlas_AD CellOracle

#################################
####### Decipher ############
#################################
sudo docker run -it -m 40g --memory-swap 44g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5
source("scripts/BCG/bcg_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/erp/erp_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/tnbc/tnbc_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/sepsis/sepsis_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/cord_pic/cord_pic_2_decipher_pipeline_v1_modularized_jan_2025.R") 
source("scripts/5yr_pic/5yr_pic_decipher_pipeline_v1_modularized_jan_2025.R") 
source("scripts/lupus/lupus_2_decipher_pipeline_v1_modularized_feb_2025.R") #done
source("scripts/covid/covid_2_decipher_pipeline_v1_modularized_feb_2025.R") 

#################################
####### Connectome ############
#################################
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-connectome@sha256:851edafd24ec10d67cd9eda8c07f611a22baf13da5667acd539bb8a390c05973
source("scripts/5yr_pic/5yr_pic_3_connectome_analysis_feb_2025.R") #done
source("scripts/cord_pic/cord_pic_3_connectome_analysis_feb_2025.R") #done
source("scripts/erp/erp_3_connectome_analysis_feb_2025.R") #requires 60 gb of RAM
source("scripts/sepsis/sepsis_3_connectome_analysis_feb_2025.R") #done
source("scripts/BCG/bcg_3_connectome_analysis_feb_2025.R") #done
source("scripts/covid/covid_3_connectome_analysis_feb_2025.R") #done
source("scripts/lupus/lupus_3_connectome_analysis_feb_2025.R") #requires 130 gb or RAM
source("scripts/tnbc/tnbc_3_connectome_analysis_feb_2025.R") #done

#################################
####### NATMI ############
#################################
#note that NATMI with docker requires downloading the databases and python scripts from the github page and placing them in a folder in the analysis directory called NATMI
#github https://github.com/asrhou/NATMI
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-natmi@sha256:1269438fa798330eba47d51ac910d76a6298fb471c6e9449685a0a1dbb2282b7 
bash scripts/5yr_pic/5yr_pic_5_natmi_analysis_feb_2025.sh
bash scripts/tnbc/tnbc_5_natmi_analysis_feb_2025.sh #done
bash scripts/sepsis/sepsis_5_natmi_analysis_feb_2025.sh #done
bash scripts/covid/covid_5_natmi_analysis_feb_2025.sh #done
bash scripts/cord_pic/cord_pic_5_natmi_analysis_feb_2025.sh #done
bash scripts/BCG/bcg_5_natmi_analysis_feb_2025.sh #done
bash scripts/erp/erp_5_natmi_analysis_feb_2025.sh #done
bash scripts/lupus/lupus_5_natmi_analysis_feb_2025.sh #running

#################################
####### NicheNet ############
#################################
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-nichenetr@sha256:146e50752019a18d4125729c3d4b09ccd88a11e188525dc48eb9dcdcb72967ec
Rscript scripts/5yr_pic/5yr_pic_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/BCG/bcg_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/cord_pic/cord_pic_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/covid/covid_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/erp/erp_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/tnbc/tnbc_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/sepsis/sepsis_4_nichenet_analysis_feb_2025.R #done
Rscript scripts/lupus/lupus_4_nichenet_analysis_feb_2025.R #done

#################################
####### LIANA+ ############
#################################

docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-liana-plus:1.0.0@sha256:d300ec7872d9a0cf8ae91fc05798f56a3aa3982657bb3883f21bcef63b8ee580 
python3 scripts/5yr_pic/5yr_pic_6_liana_plus_analysis_feb_2025.py #done
python3 scripts/BCG/bcg_6_liana_plus_analysis_feb_2025.py 
python3 scripts/cord_pic/cord_pic_6_liana_plus_analysis_feb_2025.py #done
python3 scripts/covid/covid_6_liana_plus_analysis_feb_2025.py #done
python3 scripts/erp/erp_6_liana_plus_analysis_feb_2025.py #done
python3 scripts/tnbc/tnbc_6_liana_plus_analysis_feb_2025.py #done
python3 scripts/sepsis/sepsis_6_liana_plus_analysis_feb_2025.py #done
python3 scripts/lupus/lupus_6_liana_plus_analysis_feb_2025.py #done


#################################
####### Cytosig run and process data ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-cytosig@sha256:583a450ae25f91686dbf9db9b3297d7f82f81f87b809f520daf226c3a661b11b
bash scripts/5yr_pic/5yr_pic_7_cytosig_analysis_feb_2025.sh #done
bash scripts/BCG/bcg_7_cytosig_analysis_feb_2025.sh #running
bash scripts/cord_pic/cord_pic_7_cytosig_analysis_feb_2025.sh #done
bash scripts/covid/covid_7_cytosig_analysis_feb_2025.sh #running
bash scripts/erp/erp_7_cytosig_analysis_feb_2025.sh #done
bash scripts/lupus/lupus_7_cytosig_analysis_feb_2025.sh #running
bash scripts/sepsis/sepsis_7_cytosig_analysis_feb_2025.sh #running
bash scripts/tnbc/tnbc_7_cytosig_analysis_feb_2025.sh