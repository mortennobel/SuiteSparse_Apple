#!/bin/bash

cp SuiteSparse_config/SuiteSparse_config_ios_64.mk SuiteSparse_config/SuiteSparse_config.mk
make purge
rm AMD/Lib/*.o
rm BTF/Lib/*.o
rm CAMD/Lib/*.o
rm CCOLAMD/Lib/*.o
rm CHOLMOD/Lib/*.o
rm COLAMD/Lib/*.o
rm CSparse/Lib/*.o
rm CXSparse/Lib/*.o
rm KLU/Lib/*.o
rm LDL/Lib/*.o
make library
cd SuiteSparse_config
rm -rf libsuitesparseconfig.a
make purge
make library
cd ..
cp AMD/Lib/libamd.a                               lib_ios/libamd.2.4.1_arm64.a
cp BTF/Lib/libbtf.a                               lib_ios/libbtf.1.2.1_arm64.a
cp CAMD/Lib/libcamd.a                             lib_ios/libcamd.2.4.1_arm64.a
cp CCOLAMD/Lib/libccolamd.a                       lib_ios/libcamd.2.4.1_arm64.a
cp CHOLMOD/Lib/libcholmod.a                       lib_ios/libcholmod.3.0.3_arm64.a
cp COLAMD/Lib/libcolamd.a                         lib_ios/libcolamd.2.9.1_arm64.a
cp CSparse/Lib/libcsparse.a                       lib_ios/libcsparse.3.1.4_arm64.a
cp CXSparse/Lib/libcxsparse.a                     lib_ios/libcxsparse.3.1.4_arm64.a
cp KLU/Lib/libklu.a                               lib_ios/libklu.1.3.2_arm64.a
cp LDL/Lib/libldl.a                               lib_ios/libldl.2.2.1_arm64.a
cp SuiteSparse_config/libsuitesparseconfig.a      lib_ios/libsuitesparseconfig.4.4.0_arm64.a


cp SuiteSparse_config/SuiteSparse_config_ios.mk SuiteSparse_config/SuiteSparse_config.mk
make purge
rm AMD/Lib/*.o
rm BTF/Lib/*.o
rm CAMD/Lib/*.o
rm CCOLAMD/Lib/*.o
rm CHOLMOD/Lib/*.o
rm COLAMD/Lib/*.o
rm CSparse/Lib/*.o
rm CXSparse/Lib/*.o
rm KLU/Lib/*.o
rm LDL/Lib/*.o
make library
cd SuiteSparse_config
rm -rf libsuitesparseconfig.a
make purge 
make library
cd ..
pwd
cp AMD/Lib/libamd.a                               lib_ios/libamd.2.4.1_armv7.a
cp BTF/Lib/libbtf.a                               lib_ios/libbtf.1.2.1_armv7.a
cp CAMD/Lib/libcamd.a                             lib_ios/libcamd.2.4.1_armv7.a
cp CCOLAMD/Lib/libccolamd.a                       lib_ios/libcamd.2.4.1_armv7.a
cp CHOLMOD/Lib/libcholmod.a                       lib_ios/libcholmod.3.0.3_armv7.a
cp COLAMD/Lib/libcolamd.a                         lib_ios/libcolamd.2.9.1_armv7.a
cp CSparse/Lib/libcsparse.a                       lib_ios/libcsparse.3.1.4_armv7.a
cp CXSparse/Lib/libcxsparse.a                     lib_ios/libcxsparse.3.1.4_armv7.a
cp KLU/Lib/libklu.a                               lib_ios/libklu.1.3.2_armv7.a
cp LDL/Lib/libldl.a                               lib_ios/libldl.2.2.1_armv7.a
cp SuiteSparse_config/libsuitesparseconfig.a      lib_ios/libsuitesparseconfig.4.4.0_armv7.a
cd lib_ios
./create_universal.sh
cd ..
