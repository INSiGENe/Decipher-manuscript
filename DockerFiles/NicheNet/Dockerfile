# Use r-base:4.4.1 as the base image
FROM r-base:4.4.1

# Set environment variables to prevent interactive prompts during package installations
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory inside the Docker container
WORKDIR /usr/src/app

# Install system dependencies required for R, Connectome, and nichenetr
RUN apt-get update && apt-get install -y \
    # Essential compilers and build tools
    build-essential \
    gfortran \
    libgfortran-12-dev \
    # Libraries for R packages
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
    # Install CMake for nloptr
    cmake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install necessary R packages including devtools and BiocManager
RUN R -e "install.packages(c('devtools'), repos='https://cloud.r-project.org/')"

RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')"

# Install ComplexHeatmap via BiocManager
RUN R -e "BiocManager::install('ComplexHeatmap', update=FALSE, ask=FALSE)"

# Install nichenetr from GitHub using devtools
RUN R -e "devtools::install_github('saeyslab/nichenetr', dependencies=TRUE)"
