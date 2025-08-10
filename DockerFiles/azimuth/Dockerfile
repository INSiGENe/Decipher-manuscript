FROM ebasto/decipher-manuscript-decipherc2c@sha256:7e43d263693b4c2a87a7a9459dcb1fd5ecc5a969ef84a7b3b3c2b71205efafb5

RUN R -e 'remotes::install_github("satijalab/seurat-data@v0.2.1", upgrade = "never")'

RUN R -e 'BiocManager::install("glmGamPoi",update = FALSE)'

RUN R -e 'remotes::install_github("satijalab/azimuth@release/0.4.6", upgrade = "never")'