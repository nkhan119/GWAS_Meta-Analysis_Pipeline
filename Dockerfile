# ================================================================
# GWAS / MR Pipeline — Bulletproof Build (v2.0.1)
# Author: Nadeem Khan
# Base  : Ubuntu 22.04
# ================================================================

FROM ubuntu:22.04

LABEL maintainer="Nadeem Khan"
LABEL version="2.0.1"

# ───────────────────────────────────────────────────────────────
# Environment (stability + reproducibility)
# ───────────────────────────────────────────────────────────────
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Montreal
ENV R_LIBS_SITE=/usr/local/lib/R/site-library
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
# Fixed CRAN snapshot for binary speed and stability
ENV CRAN_REPO=https://packagemanager.posit.co/cran/__linux__/jammy/2024-01-15

# ───────────────────────────────────────────────────────────────
# System dependencies
# ───────────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl git ca-certificates unzip zip pkg-config \
    libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libcairo2-dev libxt-dev libharfbuzz-dev libfribidi-dev \
    libgit2-dev libssh2-1-dev libv8-dev libudunits2-dev libgmp3-dev \
    libuv1-dev pandoc r-base r-base-dev python3 python3-pip procps \
    && rm -rf /var/lib/apt/lists/*

# ───────────────────────────────────────────────────────────────
# PLINK2
# ───────────────────────────────────────────────────────────────
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip \
    -O /tmp/plink2.zip && \
    unzip /tmp/plink2.zip -d /usr/local/bin && \
    chmod +x /usr/local/bin/plink2 && \
    rm /tmp/plink2.zip

# ───────────────────────────────────────────────────────────────
# CRAN & Bioconductor Core
# ───────────────────────────────────────────────────────────────
RUN Rscript -e " \
    options(repos = c(CRAN='${CRAN_REPO}')); \
    install.packages(c('remotes','BiocManager','R.utils','R6','rlang','cli','glue','fs','jsonlite','curl','httr','data.table','optparse','metafor','meta','qqman','CMplot','ggplot2','dplyr','magrittr'), Ncpus=parallel::detectCores()); \
    BiocManager::install(c('GenomicRanges','VariantAnnotation','AnnotationDbi'), ask=FALSE, update=FALSE); \
    "

# ───────────────────────────────────────────────────────────────
# Python deps
# ───────────────────────────────────────────────────────────────
COPY envs/requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt || true

# ───────────────────────────────────────────────────────────────
# MR ECOSYSTEM (GitHub installs)
# ───────────────────────────────────────────────────────────────
# We install dependencies IN ORDER to ensure TwoSampleMR finds them.
RUN Rscript -e " \
    options(repos = c(CRAN='${CRAN_REPO}')); \
    message('--- Installing GitHub MR Stack ---'); \
    remotes::install_github('mrcieu/ieugwasr', upgrade='never', force=TRUE); \
    remotes::install_github('rondolab/MR-PRESSO', upgrade='never', force=TRUE); \
    remotes::install_github('WSpiller/RadialMR', upgrade='never', force=TRUE); \
    remotes::install_github('gqi/MRMix', upgrade='never', force=TRUE); \
    remotes::install_github('WSpiller/MVMR', upgrade='never', force=TRUE); \
    message('--- Installing TwoSampleMR ---'); \
    remotes::install_github('MRCIEU/TwoSampleMR', upgrade='never', force=TRUE); \
    "

# ── 8. Final verification ────────────────────────────────────
# Fixed the missing braces and quotes from your previous version
RUN Rscript -e " \
    pkgs <- c('TwoSampleMR', 'MVMR', 'ieugwasr', 'MRPRESSO', 'RadialMR', 'MRMix'); \
    for (p in pkgs) { \
        if (!requireNamespace(p, quietly = TRUE)) stop(paste('Verification failed for:', p)); \
    }; \
    cat('All R packages verified successfully.\n'); \
    "

# ───────────────────────────────────────────────────────────────
# Pipeline scripts
# ───────────────────────────────────────────────────────────────
RUN mkdir -p /opt/pipeline/bin
COPY bin/ /opt/pipeline/bin/
RUN chmod +x /opt/pipeline/bin/*.py /opt/pipeline/bin/*.R || true

ENV PATH="/opt/pipeline/bin:${PATH}"
WORKDIR /data

# ───────────────────────────────────────────────────────────────
# FINAL SANITY CHECK
# ───────────────────────────────────────────────────────────────
RUN Rscript -e " \
    library(ieugwasr); \
    library(RadialMR); \
    library(TwoSampleMR); \
    cat('✔ BULLETPROOF R ENV OK\n'); \
    " && plink2 --version

CMD ["/bin/bash"]
