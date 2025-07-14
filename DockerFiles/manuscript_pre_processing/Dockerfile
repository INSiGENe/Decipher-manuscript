FROM r-base:4.4.1

# Set the working directory
WORKDIR /usr/src/app

#install system dependencies
RUN apt-get update && apt-get install -y \
    gfortran \
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
    libfontconfig1-dev \
    libhdf5-dev \
    build-essential \
    wget \
    ca-certificates \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# 2) Set renv paths *inside* container
ENV RENV_PATHS_CACHE=/renv/cache \
    RENV_PATHS_LIBRARY=/renv/library

WORKDIR /usr/src/app

# 2) Copy renv.lock, then restore packages
COPY renv.lock .

# (Optional) Install renv so that it's already available in base R.
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')" \
 && R -e "renv::restore(prompt = FALSE)"

# By default, we will start an R session in "vanilla" mode so the renv autoloader
# doesn't try to install renv behind the scenes. We’ll manage it ourselves.
# ENTRYPOINT ["R", "--vanilla"]