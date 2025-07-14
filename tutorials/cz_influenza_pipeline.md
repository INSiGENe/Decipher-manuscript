# Analysis Pipeline for the Influenza (cz_influenza) Dataset

This document outlines the full analysis pipeline for the `cz_influenza` dataset, which contains immunophenotyping data of COVID-19 and influenza.

The dataset is sourced from the CZI Science CELLxGENE database.

## 1. Dataset Download

First, download the dataset from cellxgene.

**Command:**
```bash
mkdir -p data/cz_influenza && \
wget -O data/cz_influenza/dataset.h5ad \
"https://datasets.cellxgene.cziscience.com/5f4efede-b295-4be4-aded-eb8b9a946382.h5ad"
```

## 1.1. Configuration

The following block from `config.json` is used to configure the analysis pipeline for this dataset.

```json
"cz_influenza": {
    "dataset_description":{
      "repository": "cellxgene",
      "download_url": "https://datasets.cellxgene.cziscience.com/5f4efede-b295-4be4-aded-eb8b9a946382.h5ad",
      "description": "Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19",
      "abbreviation": "cz_influenza",
      "url": "https://cellxgene.cziscience.com/collections/4f889ffc-d4bc-4748-905b-8eb9db47a2ed"
    },
    "pre_processing":{
      "input_path": "data/cz_influenza",
      "output_path": "data/cz_influenza",
      "step": "preprocessing",
      "subset_logic": {
        "disease": ["influenza", "normal"]
      },
      "min_genes_per_cell": 200
    },
    "pre_processing_for_analysis": {
      "input_path": "data/cz_influenza/preprocessing",
      "output_path": "results/cz_influenza",
      "cluster_meta_field": "cell_type",
      "condition_meta_field": "disease",
      "individual_meta_field": "donor_id",
      "case_condition": "influenza",
      "control_condition": "normal"
    },
    "cytosig": {
      "output_path": "results/cz_influenza"
    },
    "cell_oracle": {
      "data_folder": "results/cz_influenza",
      "species": "human"
    },
    "Decipher_analysis": {
      "min_cells_per_cluster_condition": 100,
      "species": "human",
      "dataset_path": "results/cz_influenza",
      "condition_name": "condition",
      "case_condition": "influenza",
      "control_condition": "normal",
      "k_parameter": 1,
      "min_meta_cells_parameter": 100
    },
    "Connectome_analysis": {
      "this_species":"human",
      "case_condition":"influenza", 
      "control_condition":"normal", 
      "dataset_path":"results/cz_influenza"
    },
    "NicheNet_analysis": {
      "case_condition": "influenza",
      "control_condition": "normal",
      "dataset_path": "results/cz_influenza",
      "this_species": "human"
    },
    "liana_plus_analysis":{
        "dataset_path": "results/cz_influenza",
        "sample_key": "donor_id",
        "groupby": "cluster",
        "condition_key": "condition"
      }
  }
```

## 2. Docker Environment Setup

The analysis pipeline relies on several Docker containers. You can either pull pre-built images from DockerHub (if available) or build them from the provided Dockerfiles.

**Build Commands:**
```bash
# For Decipher and other R-based analyses
docker build -t decipherc2c-docker:1.0.5 -f docker/Dockerfile_decipherc2c_docker_1.0.5 .

# For pre-processing steps
docker build -t manuscript_pre_processing:1.0.3 -f docker/Dockerfile_manuscript_pre_processing .
```
*Note: Other tools like CellOracle, Connectome, NicheNet, etc., use specific Docker images mentioned in their respective steps.*

## 3. Analysis Pipeline Steps

All commands should be run from the root of the project directory.

### Step 3.1: Convert AnnData to Seurat Object

This step converts the downloaded `.h5ad` file (an AnnData object) into a Seurat object, which is used in subsequent R-based analyses.

**Docker Container:** `celloracle-improved-reproducibility:latest`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest

# Run the script
python3 scripts/analysis_cellxgene_datasets/1_preprocess_h5ad.py cz_influenza
```

### Step 3.2: Generic Pre-processing

This step performs standardized pre-processing on the Seurat object. This includes filtering, normalization, and scaling, preparing the data for downstream analysis.

**Docker Container:** `manuscript_pre_processing:1.0.3`

**Command:**
```bash
# Define environment variables for renv caching
export RENV_PATHS_CACHE_HOST=/opt/local/renv/cache
export RENV_PATHS_CACHE_CONTAINER=/renv/cache

# Start the container
docker run -it --rm \
    --memory=180g --memory-swap=185g \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -v "$(pwd):/app" \
    -w /app \
    manuscript_pre_processing:1.0.3 \
    bash

# Run the script
Rscript scripts/analysis_cellxgene_datasets/2_preprocess_object_for_analysis.R cz_influenza
```

### Step 3.3: Run CytoSig

CytoSig is used to infer cell-cell communication by analyzing ligand-receptor interactions.

**Docker Container:** `data2intelligence/data2intelligence-suite`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite

# Run the script
bash scripts/analysis_cellxgene_datasets/3_cytosig_run.sh cz_influenza
```

### Step 3.4: Run CellOracle

CellOracle is used for the inference of gene regulatory networks (GRNs).

**Docker Container:** `celloracle-improved-reproducibility:latest`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest

# Run the script
python3 scripts/analysis_cellxgene_datasets/5_cell_oracle.py cz_influenza
```

### Step 3.5: Run Decipher

This is the main step of the pipeline, where Decipher is used to score and prioritize cell-cell interactions based on intracellular and intercellular events.

**Docker Container:** `decipherc2c-docker:1.0.5`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5 bash

# Run the script
Rscript scripts/analysis_cellxgene_datasets/6_decipher_pipeline_v1_modularized.R cz_influenza
```

### Step 3.6: Run Connectome

Connectome is another tool for cell-cell interaction analysis, used here for comparison.

**Docker Container:** `ebasto/connectome:latest`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/connectome:latest bash

# Run the script
Rscript scripts/analysis_cellxgene_datasets/7_connectome_analysis.R cz_influenza
```

### Step 3.7: Run NicheNet

NicheNet is used to predict ligand-target links between interacting cells.

**Docker Container:** `ebasto/nichenetr:latest`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest bash

# Run the script
Rscript scripts/analysis_cellxgene_datasets/8_nichenet_analysis.R cz_influenza
```

### Step 3.8: Run NATMI

NATMI (Network Analysis Toolkit for Multi-omics Integration) is another CCI tool used for benchmarking.

**Docker Container:** `asrhou/natmi`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace asrhou/natmi

# Run the script
bash scripts/analysis_cellxgene_datasets/9_natmi_analysis.sh cz_influenza
```

### Step 3.9: Run LIANA+

LIANA+ is a comprehensive framework that integrates multiple CCI methods.

**Docker Container:** `ebasto/liana_plus`

**Command:**
```bash
# Start the container
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus bash

# Run the script
python3 scripts/analysis_cellxgene_datasets/10_liana_plus_analysis.py cz_influenza
```
