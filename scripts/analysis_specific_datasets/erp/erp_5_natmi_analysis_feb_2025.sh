#!/bin/bash

cd NATMI
DATASET=ERP
python3 ExtractEdges.py --species human --emFile ../results/$DATASET/natmi/data/case/em.txt --annFile ../results/$DATASET/natmi/data/case/metadata.txt --interDB lrc2p --coreNum 4 --out ../results/$DATASET/natmi/data/case/
python3 ExtractEdges.py --species human --emFile ../results/$DATASET/natmi/data/control/em.txt --annFile ../results/$DATASET/natmi/data/control/metadata.txt --interDB lrc2p --coreNum 4 --out ../results/$DATASET/natmi/data/control/
python3 DiffEdges.py --refFolder ../results/$DATASET/natmi/data/control/ --targetFolder ../results/$DATASET/natmi/data/case/ --interDB lrc2p --out ../results/$DATASET/natmi/data/diff/
