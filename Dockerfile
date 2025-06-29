# Base Ubuntu image
FROM ubuntu:latest

# Install required dependencies
RUN apt-get update && apt-get install -y \
    git cmake build-essential gfortran wget \
    && rm -rf /var/lib/apt/lists/*  
    # Clean up to reduce image size

# Install OpenMPI
WORKDIR /opt
RUN wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.6.tar.gz && \
    tar -xvzf openmpi-5.0.6.tar.gz && \
    cd openmpi-5.0.6 && ./configure --prefix=/usr/local/openmpi && make -j $(nproc) all && make install

# Set OpenMPI Environment Variables
ENV PATH="/usr/local/openmpi/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/openmpi/lib:${LD_LIBRARY_PATH}"
# Set OpenMPI and PMIx environment variables
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV PMIX_MCA_pcompress_base_silence_warning=1


# Clone and Install Hypre
WORKDIR /opt
RUN git clone https://github.com/hypre-space/hypre.git && \
    cd hypre/src && ./configure --prefix=/usr/local/hypre && make -j $(nproc) all && make install

# Set Hypre Environment Variables
ENV PATH="/usr/local/hypre/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/hypre/lib:${LD_LIBRARY_PATH}"

# Clone and Build DIVEMesh
WORKDIR /opt
RUN git clone https://github.com/REEF3D/DIVEMesh.git && \
    cd DIVEMesh && make -j $(nproc)

# Clone and Build REEF3D
WORKDIR /opt
RUN git clone https://github.com/REEF3D/REEF3D.git && \
    cd REEF3D && make -j $(nproc)

# Create simulations folder and copy binaries
RUN echo "mkdir -p /simulations && cp /opt/DIVEMesh/bin/DiveMESH /simulations && cp /opt/REEF3D/bin/REEF3D /simulations" >> ~/.bashrc

# Set PATH for simulations
ENV PATH="/simulations:${PATH}"

# Set working directory for container
WORKDIR /simulations

CMD ["/bin/bash"]