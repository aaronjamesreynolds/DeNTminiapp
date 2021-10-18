#!/bin/bash

mkdir tpls; cd tpls/

#Download Third-Party libraries libraries
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz; wget 
https://github.com/hypre-space/hypre/archive/refs/tags/v2.16.0.tar.gz; wget 
https://bit.ly/mfem-4-3; wget https://github.com/hypre-
space/hypre/archive/refs/tags/v2.16.0.tar.gz;

mv mfem-4-3 mfem-4-3.tar.gz

tar -zxvf *.tar.gz

cd hypre-2.16.0/src/ 
./configure --disable-fortran 
make -j 1 
cd ../.. 
ln -s hypre-2.16.0 hypre

cd metis-4.0.3 
make 
cd .. 
ln -s metis-4.0.3 metis-4.0

module load cuda

cd mfem-4.3 
make pcuda -j 1

cd ../..
