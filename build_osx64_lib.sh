#!/bin/bash

cp SuiteSparse_config/SuiteSparse_config_Mac.mk SuiteSparse_config/SuiteSparse_config.mk

make purge
make
cp AMD/Lib/libamd.a                               lib/libamd.2.4.1_x86_64.a
cp BTF/Lib/libbtf.a                               lib/libbtf.1.2.1_x86_64.a
cp CAMD/Lib/libcamd.a                             lib/libcamd.2.4.1_x86_64.a
cp CCOLAMD/Lib/libccolamd.a                       lib/libcamd.2.4.1_x86_64.a
cp CHOLMOD/Lib/libcholmod.a                       lib/libcholmod.3.0.3_x86_64.a
cp COLAMD/Lib/libcolamd.a                         lib/libcolamd.2.9.1_x86_64.a
cp CSparse/Lib/libcsparse.a                       lib/libcsparse.3.1.4_x86_64.a
cp CXSparse/Lib/libcxsparse.a                     lib/libcxsparse.3.1.4_x86_64.a
cp KLU/Lib/libklu.a                               lib/libklu.1.3.2_x86_64.a
cp LDL/Lib/libldl.a                               lib/libldl.2.2.1_x86_64.a
cp SuiteSparse_config/libsuitesparseconfig.a      lib/libsuitesparseconfig.4.4.0_x86_64.a

cd lib
./create_universal.sh
cd ..
