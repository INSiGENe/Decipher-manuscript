# Manuscript Analysis Tutorial

This tutorial describes the computational analysis pipeline used in the manuscript. The analysis is orchestrated through a series of shell scripts that run different analysis tools, each within its own Docker container. The parameters for each analysis are controlled by a central `config.json` file.

## 1. Setup

### 1.1. Docker Images

The analysis requires several Docker images. You can either build them from the provided Dockerfiles or pull them from a Docker registry (TODO: add dockerhub links).

To build the images locally, use the following commands:

```bash
docker build -t decipherc2c-docker:1.0.5 -f Dockerfile_decipherc2c_docker_1.0.5 .
docker build -t manuscript_pre_processing:1.0.3 -f Dockerfile_manuscript_pre_processing .
# TODO: Add other docker build commands
```

### 1.2. Datasets

The datasets used in the manuscript are downloaded from public repositories. The `aws_commands/download_datasets.sh` script contains the commands to download and preprocess the data.

The main data sources are:
- cellxgene
- Gene Expression Omnibus (GEO)

To download the data, run the `download_datasets.sh` script:
```bash
bash aws_commands/download_datasets.sh
```

## 2. Pre-processing

The pre-processing steps prepare the raw data for the main analysis pipeline. These steps are detailed in the `aws_commands/analysis_commands.sh` script.

The key pre-processing steps are:
1.  **Annotation:** Datasets without annotations are annotated using Azimuth.
2.  **Format Conversion:** AnnData objects are converted to Seurat objects.
3.  **Custom Pre-processing:** Specific datasets like SevMildCOVID undergo custom pre-processing.
4.  **Generic Pre-processing:** All datasets are processed through a generic pre-processing pipeline.
5.  **scCODA Analysis:** Compositional analysis is performed using scCODA.

## 3. Analysis Pipeline

The main analysis pipeline is executed by the `aws_commands/analysis_commands.sh` script. This script iterates through the datasets defined in `manuscript_scripts/config.json` and runs a series of analysis tools.

The analysis pipeline consists of the following tools:

- Cytosig
- CellOracle
- Decipher
- Connectome
- NicheNet
- NATMI
- LIANA+

The generic command to run an analysis is:
```bash
# Example for CellOracle
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
taskset -c 0-3 python3 scripts/cell_oracle_apr_2025.py <dataset_key>
```
Where `<dataset_key>` is the name of the dataset as defined in `manuscript_scripts/config.json`.

## 4. Analysis Tools

This section will provide a more detailed description of each analysis tool and how it is used in the pipeline.

This section provides a more detailed description of each analysis tool and how it is used in the pipeline. The `5yr_pic` dataset is used as a representative example.

### 4.1. Pre-processing

The pre-processing script `5yr_pic_0_pre_processing_feb_2025_dup_2.R` reads the raw data, filters it, and prepares it for the downstream analysis. The parameters for this script are read from the `config.json` file.

**Key steps:**
- Load the Seurat object.
- Filter out cells with low gene counts.
- Subset the data based on the `disease` condition specified in the `config.json`.
- Save the processed Seurat object.

### 4.2. CellOracle

CellOracle is used to infer gene regulatory networks. The script `5yr_pic_1_cell_oracle_feb_2025_dup_2.py` runs the CellOracle analysis.

**Key steps:**
- Load the pre-processed Seurat object.
- Create a CellOracle object.
- Perform dimensionality reduction and clustering.
- Infer the gene regulatory network.
- Save the CellOracle object.

**Command:**
```bash
docker run -it -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility:latest
taskset -c 0-3 python3 manuscript_scripts/5yr_pic/5yr_pic_1_cell_oracle_feb_2025_dup_2.py 5yr_pic
```

### 4.3. Decipher

Decipher is the core tool of this analysis pipeline. The script `5yr_pic_2_decipher_pipeline_v1_modularized_feb_2025_dup_2.R` runs the Decipher analysis.

**Key steps:**
- Load the pre-processed Seurat object and the CellOracle object.
- Run the Decipher pipeline to score cell-cell interactions.
- Save the Decipher results.

**Command:**
```bash
docker run -it -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.5 bash
Rscript manuscript_scripts/5yr_pic/5yr_pic_2_decipher_pipeline_v1_modularized_feb_2025_dup_2.R 5yr_pic
```

### 4.4. Connectome

Connectome is used to identify significantly interacting cell pairs. The script `5yr_pic_3_connectome_analysis_feb_2025.R` runs the Connectome analysis.

**Command:**
```bash
docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/connectome:latest bash
Rscript manuscript_scripts/5yr_pic/5yr_pic_3_connectome_analysis_feb_2025.R 5yr_pic
```

### 4.5. NicheNet

NicheNet is used to predict ligand-receptor interactions that drive gene expression changes. The script `5yr_pic_4_nichenet_analysis_feb_2025.R` runs the NicheNet analysis.

**Command:**
```bash
docker run --rm --memory=60g --memory-swap=62g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr:latest Rscript manuscript_scripts/5yr_pic/5yr_pic_4_nichenet_analysis_feb_2025.R 5yr_pic
```

### 4.6. NATMI

NATMI is another tool for inferring cell-cell interactions. The script `5yr_pic_5_natmi_analysis_feb_2025.sh` runs the NATMI analysis.

**Command:**
```bash
docker run -it -m 30g --memory-swap 30g -v "$(pwd):/workspace" -w /workspace asrhou/natmi bash manuscript_scripts/5yr_pic/5yr_pic_5_natmi_analysis_feb_2025.sh 5yr_pic
```

### 4.7. LIANA+

LIANA+ is a tool for cell-cell communication analysis. The script `5yr_pic_6_liana_plus_analysis_feb_2025.py` runs the LIANA+ analysis.

**Command:**
```bash
docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus \
    python3 manuscript_scripts/5yr_pic/5yr_pic_6_liana_plus_analysis_feb_2025.py 5yr_pic
```

### 4.8. Cytosig

Cytosig is used to identify cell-cell interactions that are associated with a specific phenotype. The script `5yr_pic_7_cytosig_analysis_feb_2025.sh` runs the Cytosig analysis.

**Command:**
```bash
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite
bash manuscript_scripts/5yr_pic/5yr_pic_7_cytosig_analysis_feb_2025.sh 5yr_pic
```

