# ubuntu base image
FROM ubuntu:latest

# install required dependencies
RUN apt-get update
RUN apt-get -y install git cmake build-essential gfortran mpich

# clone Hypre from Github and install it
RUN git clone https://github.com/hypre-space/hypre.git
RUN cd hypre/src && ./configure --prefix=/usr/local/hypre && make all -j 8 && make install
# adding paths for hypre
RUN export PATH=/usr/local/hypre/include:/usr/local/hypre/lib:$PATH

# clone DIVEMesh and REEF3D from Github and make files
RUN git clone https://github.com/REEF3D/DIVEMesh.git
RUN cd DIVEMesh && make -j 8 && cd ..
RUN git clone https://github.com/REEF3D/REEF3D.git
RUN cd REEF3D && make -j 8

# create simulations folder and copy binaries into working folder
RUN mkdir simulations
RUN echo "export PATH=/simulations:$PATH && cp ./DIVEMesh/bin/DiveMESH ./simulations && cp ./REEF3D/bin/REEF3D ./simulations && cd simulations" >> ~/.bashrc