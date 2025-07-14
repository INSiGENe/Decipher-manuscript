# Use r-base:4.4.1 as the base image
FROM r-base:4.4.1

# Set the working directory
WORKDIR /usr/src/app

# Install system dependencies required for R and nichenetr
RUN apt-get update && apt-get install -y \
    gfortran \
    libgfortran-12-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libwebp-dev \
    pkg-config \
    git \                   
    # Required for devtools to install from GitHub
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Install necessary R packages including devtools and BiocManager
RUN R -e "install.packages(c('devtools'), repos='https://cloud.r-project.org/')"

RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')"


# Install ComplexHeatmap via BiocManager
RUN R -e "BiocManager::install('ComplexHeatmap', update=FALSE, ask=FALSE)"

# Install Connectome from GitHub using devtools
RUN R -e "devtools::install_github('msraredon/Connectome', ref = 'master', dependencies=TRUE)"





