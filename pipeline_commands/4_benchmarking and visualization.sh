docker pull ebasto/decipher-manuscript-decipherc2c@sha256:001305292bc84a0dfee0799a6bbb0f687c3d8bb2e58b2092751bbb8982a38d46

sudo docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:001305292bc84a0dfee0799a6bbb0f687c3d8bb2e58b2092751bbb8982a38d46
source("scripts/benchmarking_and_figure_scripts/1_load_all_results.r")
source("scripts/benchmarking_and_figure_scripts/figure_2.R")
source("scripts/benchmarking_and_figure_scripts/figure_3.R")
source("scripts/benchmarking_and_figure_scripts/figure_4.R")
source("scripts/benchmarking_and_figure_scripts/figure_5.R")


sudo docker run -it -v "$(pwd):/workspace" -w /workspace ebasto/decipher-manuscript-decipherc2c@sha256:001305292bc84a0dfee0799a6bbb0f687c3d8bb2e58b2092751bbb8982a38d46
source("scripts/benchmarking_and_figure_scripts/2a_robustness_benchmarking.R")
