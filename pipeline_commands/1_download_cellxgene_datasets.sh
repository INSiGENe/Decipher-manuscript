#navigate to analysis directory

#################################
####### Download datasets ############
#################################

#Run these commands at the root of the analysis folder, datasets will be downloaded to their own data folder.

######################################
# CellXGene
######################################

#human kidney
mkdir -p data/cz_human_kidney_v1.5 && \
wget -O data/cz_human_kidney_v1.5/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/f5b6d620-76df-45c5-9524-e5631be0e44a.h5ad"

#cz-influenza
mkdir -p data/cz_influenza && \
wget -O data/cz_influenza/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/5f4efede-b295-4be4-aded-eb8b9a946382.h5ad"

#hpap islets
mkdir -p data/cz_hpap_t1d_islets && \
wget -O data/cz_hpap_t1d_islets/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/f89a618b-fe4b-404e-bd39-7c574529b1f5.h5ad"

#chrons
mkdir -p data/cz_dev_gut_crohns && \
wget -O data/cz_dev_gut_crohns/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/518edb4f-1e86-4af3-99e4-0a599a94e9fe.h5ad"

#hnscc hpv
mkdir -p data/cz_hnscc_hpv && \
wget -O data/cz_hnscc_hpv/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/85981f97-aaa8-4c45-b920-7522727ea58c.h5ad"

#bronchial
mkdir -p data/cz_cf_bronchial_biopsy && \
wget -O data/cz_cf_bronchial_biopsy/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/a6657cae-5daa-45cd-b1b4-cf08a07d3a7e.h5ad"

#### --------------------------------- ####
####   Custom download (SevMildCOVID)  ####
#### --------------------------------- ####

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
      echo "$filename not valid. Removing."
      rm "$outpath"
    fi
  done
done


# Shared Files to download
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
####### Legacy datasets (still used but more custom pre-processing) #########
#################################

#what about other ones?
mkdir -p "data/erp"
mkdir -p "data/tnbc"
#this data sits behind a login so needs to be downloaded manually
#add 1863-counts_cells_cohort1.rds and 1872-BIOKEY_metaData_cohort1_web.csv.xls to the folder above

#PIC stimulation
mkdir -p "data/pic_5yr"
mkdir -p "data/pic_cord"
#for both these datasets, add the Zenodo-stored seurat_object_integrated_annotated.rds file (which corresponds to data from Read et al.) to the folders above

# BCG baccination
mkdir -p "data/bcg"
#place the GEO-downloaded GSE244126_RAW (unzipped folder) in the folder above

#covid (pfizer vaccination) 
mkdir -p "data/covid"
#add GSE171964_geo_pheno_v2.csv.gz unzipped and with name GSE171964_geo_pheno_v2.csv to the folder above
mkdir -p "data/covid/GSE171964_unzipped"
#add the GEO-downloaded GSE171964_barcodes_v2.tsv.gz, GSE171964_countsmatrix_v2.mtx.gz and GSE171964_feats_v2.tsv.gz (unzipped and renamed to matrix.mtx, barcodes.tsv, features.tsv, respectively to the folder above)

mkdir -p "data/lupus"
# place the cellxgene .h5ad file in the folder above and rename to lupus_sc_data.h5ad

#sepsis
mkdir -p "data/sepsis" && \
wget --show-progress -O "data/sepsis/matrix.csv.gz" \
  "https://singlecell.broadinstitute.org/single_cell/data/public/SCP548/an-immune-cell-signature-of-bacterial-sepsis-patient-pbmcs?filename=scp_gex_matrix_raw.csv.gz" && \
wget --show-progress -O "data/sepsis/meta_data.txt" \
  "https://singlecell.broadinstitute.org/single_cell/data/public/SCP548/an-immune-cell-signature-of-bacterial-sepsis-patient-pbmcs?filename=scp_meta_updated.txt"



##Not done, but for future reference, can add lupus here
#mkdir -p data/lupus && \
#wget -O data/lupus/dataset.h5ad \
#"https://datasets.cellxgene.cziscience.com/6e65a13e-2a45-4340-9442-bd2c41a01f17.h5ad"
