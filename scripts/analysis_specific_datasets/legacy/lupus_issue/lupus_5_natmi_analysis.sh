#!/bin/bash

cd NATMI
DATASET=lupus
python3 ExtractEdges.py --species human --emFile ../manuscript_analysis/$DATASET/natmi/data/case/em.txt --annFile ../manuscript_analysis/$DATASET/natmi/data/case/metadata.txt --interDB lrc2p --coreNum 4 --out ../manuscript_analysis/$DATASET/natmi/data/case/
python3 ExtractEdges.py --species human --emFile ../manuscript_analysis/$DATASET/natmi/data/control/em.txt --annFile ../manuscript_analysis/$DATASET/natmi/data/control/metadata.txt --interDB lrc2p --coreNum 4 --out ../manuscript_analysis/$DATASET/natmi/data/control/
python3 DiffEdges.py --refFolder ../manuscript_analysis/$DATASET/natmi/data/control/ --targetFolder ../manuscript_analysis/$DATASET/natmi/data/case/ --interDB lrc2p --out ../manuscript_analysis/$DATASET/natmi/data/diff/
