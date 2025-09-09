docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-celloracle@sha256:419f15c53249e4c09a0c361b8e9ab0857d46c1e797d0f3c40af35a1ce583c1b1

#generic analysis command
python3 scripts/sample_analysis/cell_oracle.py