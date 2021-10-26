nvcc -O3 -std=c++11 -ccbin mpicxx -I../tpls/mfem-4.2 -I../tpls/hypre/src/hypre/include miniapp.cpp -L../tpls/mfem-4.2 -lmfem -L../tpls/hypre/src/hypre/lib -lHYPRE -lcusparse -lcurand -L../tpls/metis-4.0 -lmetis -lcusparse -lrt -o DeNTapp-gpu

