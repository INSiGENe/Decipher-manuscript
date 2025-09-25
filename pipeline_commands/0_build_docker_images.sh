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

#
docker run --rm -it -v "$(pwd):/home/project" -w /home/project ebasto/decipherc2c@sha256:9ec8bbb7eca692856e9fd7c50346c78eac6ae24da8bd2e87b1f36998465bd1b0 bash


