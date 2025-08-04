docker pull ebasto/decipher-manuscript-decipherc2c@sha256:001305292bc84a0dfee0799a6bbb0f687c3d8bb2e58b2092751bbb8982a38d46

sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:001305292bc84a0dfee0799a6bbb0f687c3d8bb2e58b2092751bbb8982a38d46
source("scripts/benchmarking_and_figure_scripts/1_load_all_results.r")


sudo docker run -it -m 20g --memory-swap 24g -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:001305292bc84a0dfee0799a6bbb0f687c3d8bb2e58b2092751bbb8982a38d46
source("scripts/benchmarking_and_figure_scripts/2a_robustness_benchmarking.R")