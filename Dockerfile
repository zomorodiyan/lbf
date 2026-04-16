# =============================================================================
# Dockerfile for laserbeamFoam (lbf3)
# Based on OpenFOAM v2506 (opencfd/openfoam-dev)
# =============================================================================

FROM opencfd/openfoam-dev:2506

USER root

# Install build dependencies (including VTK for LIGGGHTS)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        gfortran \
        libfftw3-dev \
        libjpeg-dev \
        libpng-dev \
        libvtk9-dev \
        libopenmpi-dev \
        openmpi-bin \
    && rm -rf /var/lib/apt/lists/*

# ---- Build LIGGGHTS (optional – comment out this block if not needed) ------
RUN git clone --depth 1 https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git \
        /opt/LIGGGHTS-PUBLIC \
    && cd /opt/LIGGGHTS-PUBLIC/src \
    && make -j"$(nproc)" auto \
    && ln -s /opt/LIGGGHTS-PUBLIC/src/lmp_auto /usr/local/bin/liggghts

# ---- Build laserbeamFoam ---------------------------------------------------
COPY . /opt/laserbeamFoam

RUN bash -c '\
    source /usr/lib/openfoam/openfoam2506/etc/bashrc \
    && cd /opt/laserbeamFoam \
    && ./Allwmake -j \
    '

# Ensure the OpenFOAM environment is sourced for every bash session
RUN echo 'source /usr/lib/openfoam/openfoam2506/etc/bashrc' >> /etc/bash.bashrc

# Put compiled solver/utility binaries on PATH
ENV PATH="/root/OpenFOAM/user-v2506/platforms/linux64GccDPInt32Opt/bin:${PATH}"
ENV LD_LIBRARY_PATH="/root/OpenFOAM/user-v2506/platforms/linux64GccDPInt32Opt/lib"

# Default working directory: tutorials
WORKDIR /opt/laserbeamFoam/tutorials

# Default command: interactive shell with OpenFOAM sourced
CMD ["bash", "-l"]
