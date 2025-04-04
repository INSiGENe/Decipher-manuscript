cd projects/Manuscript_jan_2025
#################################
####### Build Docker images ############
#################################
docker build -t decipherc2c-docker:1.0.5 -f Dockerfile_decipherc2c_docker_1.0.5 .
docker build -t manuscript_pre_processing:1.0.3 -f Dockerfile_manuscript_pre_processing .

#################################
####### Download datasets ############
#################################
mkdir -p data/cz_placenta_infection && \
wget -O data/cz_placenta_infection/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/50bbe1a2-5f27-47f5-a809-046459a4ae5e.h5ad"


mkdir -p data/human_kidney_v1.5 && \
wget -O data/human_kidney_v1.5/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/f5b6d620-76df-45c5-9524-e5631be0e44a.h5ad"

mkdir -p data/cz_periheart && \
wget -O data/cz_periheart/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/66c82e1b-e3ce-48dc-b1db-52546dbd4e44.h5ad"

mkdir -p data/cz_carebank && \
wget -O data/cz_carebank/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/e67f1a92-0371-4657-b15e-a4934f9ab733.h5ad"

mkdir -p data/cz_rcc && \
wget -O data/cz_rcc/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/bf28d870-0750-443d-bb54-ec664b8f40c3.h5ad"

mkdir -p data/cz_rcc && \
wget -O data/cz_rcc/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/bf28d870-0750-443d-bb54-ec664b8f40c3.h5ad"

mkdir -p data/cz_influenza && \
wget -O data/cz_influenza/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/5f4efede-b295-4be4-aded-eb8b9a946382.h5ad"

mkdir -p data/cz_afib_macrophages && \
wget -O data/cz_afib_macrophages/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/bca7945f-69e7-4bab-92d3-90e6af99c7ac.h5ad"

mkdir -p data/cz_hpap_t1d_islets && \
wget -O data/cz_hpap_t1d_islets/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/f89a618b-fe4b-404e-bd39-7c574529b1f5.h5ad"

mkdir -p data/cz_dev_gut_crohns && \
wget -O data/cz_dev_gut_crohns/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/518edb4f-1e86-4af3-99e4-0a599a94e9fe.h5ad"

mkdir -p data/cz_hnscc_hpv && \
wget -O data/cz_hnscc_hpv/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/85981f97-aaa8-4c45-b920-7522727ea58c.h5ad"

mkdir -p data/cz_ra_pbmc && \
wget -O data/cz_ra_pbmc/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/ac9c13da-7134-4d09-8086-d0933cbdba41.h5ad"

mkdir -p data/cz_cf_bronchial_biopsy && \
wget -O data/cz_cf_bronchial_biopsy/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/a6657cae-5daa-45cd-b1b4-cf08a07d3a7e.h5ad"


#severe-mild covid
#!/bin/bash

# Parent directory name
parent_dir="data/SevMilCOVID"
mkdir -p "$parent_dir"

# Base URL
base_url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155673&format=file&file="

# Only real RNA samples available in GEO
samples=(cov01 cov02 cov03 cov04 cov07 cov08 cov09 cov10 cov11 cov12 cov17 cov18)

# Loop from cov01 to cov18
for sample_id in "${samples[@]}"; do
  folder="${parent_dir}/${sample_id}"
  mkdir -p "$folder"

  for suffix in "barcodes.tsv.gz" "matrix.mtx.gz"; do
    filename="GSE155673_${sample_id}_${suffix}"
    url="${base_url}${filename}"
    outpath="${folder}/${filename}"

    echo "Downloading $filename into $folder"
    curl -s -L "$url" -o "$outpath"

    # Check if it downloaded a valid gzip
    if file "$outpath" | grep -qv "gzip compressed data"; then
      echo "❌ $filename not valid. Removing."
      rm "$outpath"
    fi
  done
done


# common Files to download
files=(
  "GSE155673_barcodes.tsv.gz"
  "GSE155673_features.tsv.gz"
  "GSE155673_totalseq_a_hto.csv.gz"
)

# Download each file into the parent directory
for filename in "${files[@]}"; do
  url="${base_url}${filename}"
  echo "Downloading $filename into $parent_dir"
  curl -s -L "$url" -o "${parent_dir}/${filename}"
done

cd data/SevMildCOVID

# Unzip all .gz files in the root directory and delete the .gz files
echo "Unzipping root-level files..."
gunzip *.gz

# Loop through each subfolder starting with 'cov' and unzip files inside, deleting .gz files
for d in cov*/; do
  echo "Unzipping files in $d..."
  gunzip "${d}"*.gz
done


#################################
####### Convert anndata objects to R-based objects (Seurat) ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
python3 scripts/preprocess_ann_object.py dataset_key
#example
python3 scripts/preprocess_ann_object.py cz_cf_bronchial_biopsy

#################################
####### Process Initial Seurat object for downstream analyses ############
#################################
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

#run custom pre-processing if exists first, this sets it up for preprocess_object_for_analysis
#SevMilCovid (Severe and Moderate Covid-19 datasets)
Rscript scripts/custom_pre_processing_SevMilCovid.r SevCOVID
Rscript scripts/custom_pre_processing_SevMilCovid.r MilCOVID

#then run
Rscript scripts/preprocess_object_for_analysis.R dataset_key

#################################
####### Move results to analysis folder ############
#################################
sudo mv pre_processing_test/results/dataset_key Manuscript_jan_2025/results/

#example
sudo mv pre_processing_test/results/MilCOVID Manuscript_jan_2025/results/

#################################
####### Run cytosig analysis ############
#################################
cd projects/Manuscript_jan_2025
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite
bash scripts/cytosig_run.sh dataset_key

#example
bash scripts/cytosig_run.sh SevCOVID

#################################
####### Run CellOracle analysis ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
taskset -c 0-3 python3 scripts/cell_oracle_apr_2025.py dataset_key

#example
taskset -c 0-1 python3 scripts/cell_oracle_apr_2025.py SevCOVID
taskset -c 2-3 python3 scripts/cell_oracle_apr_2025.py MilCOVID


#################################
####### Decipher ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5 bash
Rscript scripts/decipher_pipeline_v1_modularized.R dataset_key

#example
Rscript scripts/decipher_pipeline_v1_modularized.R cz_dev_gut_crohns


#################################
####### Connectome ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/connectome:latest bash
Rscript scripts/connectome_analysis.R MilCOVID

#################################
####### NicheNet ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest bash
Rscript scripts/nichenet_analysis.R MilCOVID

#################################
####### NATMI ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace asrhou/natmi  
bash scripts/natmi_analysis.sh cz_cf_bronchial_biopsy

#################################
############# LIANA+ ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus bash
python3 scripts/liana_plus_analysis.py cz_cf_bronchial_biopsy




























#################################
####### Run Decipher analysis ############
#################################
sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5

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
taskset -c 0-3 python cell_oracle_mar_2025.py SkinAtlas_AD CellOracle
docker run -it -m 60g --memory-swap 62g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest

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

#NicheNet
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/5yr_pic/5yr_pic_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/BCG/bcg_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/cord_pic/cord_pic_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/covid/covid_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/erp/erp_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/tnbc/tnbc_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/sepsis/sepsis_4_nichenet_analysis_feb_2025.R #done
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript scripts/lupus/lupus_4_nichenet_analysis_feb_2025.R #done

docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/5yr_pic/5yr_pic_6_liana_plus_analysis_feb_2025.py #done
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/BCG/bcg_6_liana_plus_analysis_feb_2025.py 
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/cord_pic/cord_pic_6_liana_plus_analysis_feb_2025.py #done
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/covid/covid_6_liana_plus_analysis_feb_2025.py #done
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/erp/erp_6_liana_plus_analysis_feb_2025.py #done
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/tnbc/tnbc_6_liana_plus_analysis_feb_2025.py #done
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/sepsis/sepsis_6_liana_plus_analysis_feb_2025.py #done
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 scripts/lupus/lupus_6_liana_plus_analysis_feb_2025.py #done


#################################
####### Cytosig run and process data ############
#################################
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite
bash scripts/5yr_pic/5yr_pic_7_cytosig_analysis_feb_2025.sh #done
bash scripts/BCG/bcg_7_cytosig_analysis_feb_2025.sh #running
bash scripts/cord_pic/cord_pic_7_cytosig_analysis_feb_2025.sh #done
bash scripts/covid/covid_7_cytosig_analysis_feb_2025.sh #running
bash scripts/erp/erp_7_cytosig_analysis_feb_2025.sh #done
bash scripts/lupus/lupus_7_cytosig_analysis_feb_2025.sh #running
bash scripts/sepsis/sepsis_7_cytosig_analysis_feb_2025.sh #running
bash scripts/tnbc/tnbc_7_cytosig_analysis_feb_2025.sh



docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite
bash scripts/cytosig_run.sh cz_placenta_infection #done
bash scripts/cytosig_run.sh cz_rcc #running
bash scripts/cytosig_run.sh cz_human_kidney_v1.5  #running


python3 scripts/cytosig_run.sh cz_influenza
python3 scripts/cytosig_clean.py cz_influenza