# Decipher
A computational pipeline to extract context-specific mechanistic insights from single-cell profiles

# Decipher Pipeline

This repository contains the scripts and configuration for running the **Decipher** analysis pipeline.  
The private release of the scripts and reference data required to reproduce the analysis in the Decipher manuscript is available via **Zenodo**.

## Getting Started

### 1. Download the Repository and Data
If reproducing the results from the Decipher manuscript, please download the version from **Zenodo**.  

---

## Pipeline Overview

The pipeline is orchestrated using pipeline commands in the **pipeline_commands** folder, located in the root of the repository. These are numbered `1` to `4`:

1. **`download-datasets.sh`**  
   - Downloads datasets from CELLxGENE and GEO (e.g., severe vs mild COVID-19) and set up directory structures for manual-analysis/download (MAD) datasets.  
   - MAD datasets download instructions are documented in the paper’s **Methods** and **Data Availability** sections.

2. **`dynamic_commands_analysis`**  
   - Uses `scripts/config.json` to pull dataset-specific parameters dynamically.  
   - Pulls Docker images for preprocessing and multiple benchmarking methods:  
     - **Decipher** (R-based)
     - **Cytosig** (R-based)
     - **CellOracle** (python-based, prerequisite for Decipher)
     - **Connectome**
     - **NATMI**
     - **LIANA+**  
     - **NicheNet**
     - **Cytosig**
   - Preprocessing steps convert Python `.h5ad` files to R-compatible formats.

3. **`static_commands_analysis`**  
   - For datasets without `config.json` parameters, commands are defined directly in the scripts, hence each dataset and step has a unique script. 
   - Each section begins with a Docker pull command, followed by a set of independent commands that can be run sequentially or in parallel.

4. **figures**  
   - Generates figures for the manuscript.  
   - Each figure set is grouped under a Docker pull command:
     - **Figure 2** results must be run sequentially.
     - **Figures 3–5** and other figures can be run independently.

---

## Supporting Files and Directories

- **`reference_data/`**  
  Contains required reference data for analysis. Large files (e.g., NicheNet reference data) are excluded from Git tracking and so are not included in the Github repository, but are available via Zendo

- **`scripts/`**  
  Holds individual analysis scripts for preprocessing, benchmarking, and figure generation.

- **`scripts/config.json`**  
  Defines parameters for datasets processed dynamically.  
  Structure:
  - Top-level keys: dataset names
  - Nested keys: preprocessing and method-specific settings (e.g., Cytosig, Loracle, Cipher)

- **`dockerfiles`**
  Contains select dockerfile commands to rebuild docker images

- **`sample_analysis/`**  
  Contains sample analysis data (to be reviewed for updates).

- **`tutorials/`** *(need to be updated)*

- **`NATMI/`**  
  Python functions and databases used for running NATMI analysis 

---

## Running the Pipeline

Please step through each bash file within pipeline_commands.

---

