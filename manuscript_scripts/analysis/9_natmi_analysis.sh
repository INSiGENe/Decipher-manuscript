#!/bin/bash
# Usage: ./run_cytosig.sh cz_rcc
set -e

# Dataset key passed as argument
dataset_key=$1

cd NATMI

case_out="../results/${dataset_key}/natmi/data/case/Edges_lrc2p.csv"
control_out="../results/${dataset_key}/natmi/data/control/Edges_lrc2p.csv"

# Run ExtractEdges for case
if [ ! -f "$case_out" ]; then
  echo "Running ExtractEdges for case..."
  python3 ExtractEdges.py \
    --species human \
    --emFile ../results/${dataset_key}/natmi/data/case/em.txt \
    --annFile ../results/${dataset_key}/natmi/data/case/metadata.txt \
    --interDB lrc2p \
    --coreNum 4 \
    --out ../results/${dataset_key}/natmi/data/case/
else
  echo "Case output already exists, skipping."
fi

# Run ExtractEdges for control
if [ ! -f "$control_out" ]; then
  echo "Running ExtractEdges for control..."
  python3 ExtractEdges.py \
    --species human \
    --emFile ../results/${dataset_key}/natmi/data/control/em.txt \
    --annFile ../results/${dataset_key}/natmi/data/control/metadata.txt \
    --interDB lrc2p \
    --coreNum 4 \
    --out ../results/${dataset_key}/natmi/data/control/
else
  echo "Control output already exists, skipping."
fi

# Run DiffEdges
python3 DiffEdges.py --refFolder ../results/${dataset_key}/natmi/data/control/ --targetFolder ../results/${dataset_key}/natmi/data/case/ --interDB lrc2p --out ../results/${dataset_key}/natmi/data/diff/
