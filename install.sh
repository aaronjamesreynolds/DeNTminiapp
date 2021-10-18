#!/bin/bash
module load cmake
module load cuda
module load hypre
module load metis

mkdir build
cd build


#Download necessary libraries
wget -O mfem-4.3.tgz "https://bit.ly/mfem-4-3"
wget -O hypre-2.23.0.tar.gz "https://github.com/hypre-space/hypre/archive/refs/tags/v2.23.0.tar.gz"
wget -O metis-5.1.0.tar.gz "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"
wget -O glvis-4.1.tgz "https://bit.ly/glvis-4-1"


#build mfem serial
tar -zxvf mfem-4.3.tgz
cd mfem-4.3
make serial -j 8
cd ..

##build hypre
tar -zxvf hypre-2.23.0.tar.gz
cd hypre-2.23.0/src/
./configure --disable-fortran
make install -j
cd ../..
ln -s hypre-2.23.0 hypre


##Build METIS 5
tar zvxf metis-5.1.0.tar.gz
cd metis-5.1.0
make config ; make
mkdir lib
ln -s ../build/Linux-x86_64/libmetis/libmetis.a lib
cd ..

#build mfem
cd mfem-4.3
#make pcuda -j 8 MFEM_USE_METIS_5=YES METIS_DIR=../metis-5.1.0
#make cuda -j 8
cd ..

#build glvis
tar -zxvf glvis-4.1.tgz
cd glvis-4.1
make MFEM_DIR=../mfem-4.3 -j
cd ..
