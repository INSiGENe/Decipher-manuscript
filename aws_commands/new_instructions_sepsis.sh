#activate celloracle env in a net terminal window and run
#change to analysis directory
python3 scripts/sepsis_1_cell_oracle.py

#activate natmi env
bash scripts/sepsis_5_natmi_analysis.sh

#activate liana env
python3 scripts/sepsis_6_liana_plus_analysis.py



#for cytosig
docker run -it -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite
bash scripts/sepsis_7_cytosig_analysis.sh
python3 scripts/silent/sepsis_cytosig_process_files.py


docker run -it -v "$(pwd):/workspace" -w /workspace kenjikamimoto126/celloracle_ubuntu:0.18.0

docker run -it -m 12g --memory-swap 16g -v "$(pwd):/workspace" -w /workspace kenjikamimoto126/celloracle_ubuntu:0.18.0

