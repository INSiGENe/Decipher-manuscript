docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1

#generic analysis command
#NOTE: if you want to run multiple CellOracle analyses in parallel, constraining the number of cores per analysis prevents overflow
taskset -c 0-3 python3 scripts/analysis_cellxgene_datasets/5_cell_oracle.py dataset_key