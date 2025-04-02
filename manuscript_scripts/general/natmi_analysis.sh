#!/bin/bash
# Usage: ./run_cytosig.sh cz_rcc
set -e

# Dataset key passed as argument
dataset_key=$1

cd NATMI

python3 ExtractEdges.py --species human --emFile ../results/${dataset_key}/natmi/data/case/em.txt --annFile ../results/${dataset_key}/natmi/data/case/metadata.txt --interDB lrc2p --coreNum 4 --out ../results/${dataset_key}/natmi/data/case/
python3 ExtractEdges.py --species human --emFile ../results/${dataset_key}/natmi/data/control/em.txt --annFile ../results/${dataset_key}/natmi/data/control/metadata.txt --interDB lrc2p --coreNum 4 --out ../results/${dataset_key}/natmi/data/control/
python3 DiffEdges.py --refFolder ../results/${dataset_key}/natmi/data/control/ --targetFolder ../results/${dataset_key}/natmi/data/case/ --interDB lrc2p --out ../results/${dataset_key}/natmi/data/diff/
