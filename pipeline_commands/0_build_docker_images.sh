#build the development image
docker build -t ebasto/decipherc2c:base-dev .
#tag - sha256:30b8e41cfe0ac5d648d6a695613754b85414912ed2dec08f64bb27e9c37a6a91 

#use like this - basically uses the renv cache
export RENV_PATHS_CACHE_HOST=/opt/local/renv/cache
# The path *inside* the container that we will mount it to.
export RENV_PATHS_CACHE_CONTAINER=/renv/cache
docker run -it --rm \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -v "$(pwd):/app" \
    -w /app \
    ebasto/decipherc2c:base-dev \
    bash

#ok so we need to transfer over our R functions primarily and other scripts that we use to the 

 lockfile <- renv::lockfile_read("renv.lock")
 packages <- names(lockfile$Packages)
 renv::install(packages,exclude = c("DelayedArray"))

 renv::install("bioconductor::DelayedArray")

