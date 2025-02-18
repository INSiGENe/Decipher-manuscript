
sudo docker run -it -m 100g --memory-swap 110g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker  
sudo docker run -it -m 45g --memory-swap 50g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility



sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker:1.0.3




#run this at the head of the folder structure
sudo docker run -it -m 50g --memory-swap 55g -v "$(pwd):/workspace" -w /workspace ebasto/celloracle_mod:1.0.0
python3 analysis/scripts/aws_instance_celloracle.py


sudo docker run -it -m 8g --memory-swap 10g -v "$(pwd):/workspace" -w /workspace celloracle-improved-reproducibility

sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace ebasto/connectome
sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace ebasto/nichenetr
source("scripts/local_dev_comparison_1/4_nichenet_analysis.R")

sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace asrhou/natmi

sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace ebasto/liana_plus


sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace data2intelligence/data2intelligence-suite

sudo docker run -it -m 10g --memory-swap 12g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker  

sudo docker run -it -m 220g --memory-swap 225g -v "$(pwd):/workspace" -w /workspace decipherc2c-docker  



