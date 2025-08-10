FROM satijalab/azimuth

# Install system libraries required for devtools-related packages
RUN apt-get update && apt-get install -y \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Fix for systemfonts: preinstall cpp11 (no other packages)
RUN Rscript -e "install.packages('cpp11', repos='https://cloud.r-project.org')"


# Install devtools and remotes
RUN Rscript -e "install.packages(c('devtools', 'remotes'), repos='https://cloud.r-project.org')"

# Install seurat-data from GitHub via remotes
RUN Rscript -e "remotes::install_github('satijalab/seurat-data')" 
