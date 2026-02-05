# Use the official RStudio image as a base
# download pre-build rocker image that comes with R studio server installed. 
FROM rocker/rstudio:latest  

# Install system dependencies for python
# the -y optios gives yes automatically to any response, otherwise it defaults to no and installation fails.
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('tidyverse','caret','devtools','shiny','reticulate','pryr','umap','doParallel',
'randomForest','MLmetrics','kernlab','naivebayes','mclust'))"

# install python packages
RUN pip install --no-cache-dir numpy pandas scikit-learn dask dask[dataframe]

# Set up a working directory
WORKDIR /home/rstudio

# Set environment variables
ENV PASSWORD=1234

# Expose the RStudio port
EXPOSE 8787

# Set the entrypoint for RStudio
# /init is a initiation script included in the rocker image
CMD ["/init"]



