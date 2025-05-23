#!/bin/bash
set -eu


## Constants
#
VER="6.2"

PREFIX="${ISSM_DIR}/externalpackages/dakota/install" # Set to location where external package should be installed

# Find libgfortran and libgcc so we do not have to hardcode them
#
# TODO:
# - Move this to etc/environment.sh
# - Test if -static-libgfortran flag will avoid all of this.
# - Otherwise, refactor this to work with other gfortran installations.
#
echo "Finding libgfortran..."
LIBGFORTRAN=$(find /usr /opt -name libgfortran* 2>/dev/null | egrep -n libgfortran.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
LIBGFORTRAN_ROOT=${LIBGFORTRAN%/*}
LIBGCC=$(find ${LIBGFORTRAN_ROOT} -name libgcc* 2>/dev/null | egrep -n libgcc.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)

## Environment
#
export BLAS_LIBS="-L${BLAS_ROOT}/lib -lfblas ${LIBGFORTRAN_ROOT}/libgfortran.a ${LIBGFORTRAN_ROOT}/libquadmath.a ${LIBGCC}" # Need to export BLAS_LIBS *and* pass it as an option to CMake to ensure that external packages also find it
export DAK_BUILD=${ISSM_DIR}/externalpackages/dakota/build # DO NOT CHANGE THIS
export DAK_INSTALL=${PREFIX} # DO NOT CHANGE THIS
export DAK_SRC=${ISSM_DIR}/externalpackages/dakota/src # DO NOT CHANGE THIS
export LAPACK_LIBS="-L${LAPACK_ROOT}/lib -lflapack ${LIBGFORTRAN_ROOT}/libgfortran.a ${LIBGFORTRAN_ROOT}/libquadmath.a ${LIBGCC}" # Need to export LAPACK_LIBS *and* pass it as an option to CMake to ensure that external packages also find it
export LDFLAGS="-framework CoreFoundation"

# Cleanup
rm -rf ${DAK_BUILD} ${DAK_INSTALL} ${DAK_SRC}
mkdir -p ${DAK_BUILD} ${DAK_INSTALL} ${DAK_SRC}

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/dakota-${VER}-public.src.tar.gz" "dakota-${VER}-public-src.tar.gz"

# Unpack source
tar -zxvf dakota-${VER}-public-src.tar.gz

# Move source to 'src' directory
mv dakota-${VER}.0.src/* ${DAK_SRC}
rm -rf dakota-${VER}.0.src

# Copy customized source and configuration files to 'src' directory
cp configs/${VER}/packages/DDACE/src/Analyzer/MainEffectsExcelOutput.cpp ${DAK_SRC}/packages/DDACE/src/Analyzer
cp configs/${VER}/packages/surfpack/src/surfaces/nkm/NKM_KrigingModel.cpp ${DAK_SRC}/packages/surfpack/src/surfaces/nkm
cp configs/${VER}/packages/VPISparseGrid/src/sandia_rules.cpp ${DAK_SRC}/packages/VPISparseGrid/src
cp configs/${VER}/src/DakotaInterface.cpp ${DAK_SRC}/src
cp configs/${VER}/src/NonDLocalReliability.cpp ${DAK_SRC}/src
cp configs/${VER}/src/NonDSampling.cpp ${DAK_SRC}/src
cp configs/${VER}/static/cmake/BuildDakotaCustom.cmake ${DAK_SRC}/cmake
cp configs/${VER}/cmake/DakotaDev.cmake ${DAK_SRC}/cmake

# Copy customized source and configuration files specific to Mac to 'src' directory
cp configs/${VER}/mac/cmake/InstallDarwinDylibs.cmake ${DAK_SRC}/cmake

# Patch source
patch ${DAK_SRC}/src/dakota_data_io.hpp configs/${VER}/src/dakota_data_io.hpp.patch

# Uncomment to copy the following customized source files if using C99 or later 
# standard. If uncommented, adding -Wno-error=implicit-function-declaration 
# option to CFLAGS is not needed.
# cp configs/${VER}/mac/static/packages/DDACE/include/xtndispatch.h ${DAK_SRC}/packages/DDACE/include
# cp configs/${VER}/mac/static/packages/DDACE/src/Bose/Boselink.c ${DAK_SRC}/packages/DDACE/src/Bose
# cp configs/${VER}/mac/static/packages/DDACE/src/Bose/construct.c ${DAK_SRC}/packages/DDACE/src/Bose
# cp configs/${VER}/mac/static/packages/DDACE/src/Bose/galois.c ${DAK_SRC}/packages/DDACE/src/Bose
# cp configs/${VER}/mac/static/packages/DDACE/src/Bose/gfields.c ${DAK_SRC}/packages/DDACE/src/Bose
# cp configs/${VER}/mac/static/packages/DDACE/src/Bose/oa.c ${DAK_SRC}/packages/DDACE/src/Bose
# cp configs/${VER}/mac/static/packages/nidr/nidr.c ${DAK_SRC}/packages/nidr
# cp configs/${VER}/mac/static/packages/nidr/nidr-scanner.c ${DAK_SRC}/packages/nidr

# Disable requirement of Python 2 for TriBITS
sed -i'' -e 's|SET(PythonInterp_FIND_VERSION|#SET(PythonInterp_FIND_VERSION|' ${DAK_SRC}/packages/teuchos/cmake/tribits/package_arch/TribitsFindPythonInterp.cmake

# Configure
#
# NOTE:
# - The -w option has been added to CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, and 
#	CMAKE_Fortran_FLAGS. This should be removed for more recent versions of 
#	Dakota.
#
cd ${DAK_BUILD}
cmake \
	-DBUILD_SHARED_LIBS=OFF \
	-DBUILD_STATIC_LIBS=ON \
	-DCMAKE_C_COMPILER=${MPI_HOME}/bin/mpicc \
	-DCMAKE_C_FLAGS="-fPIC -w -Wno-error=implicit-int -Wno-error=implicit-function-declaration" \
	-DCMAKE_CXX_COMPILER=${MPI_HOME}/bin/mpicxx \
	-DCMAKE_CXX_FLAGS="-fPIC -fdelayed-template-parsing -w" \
	-DCMAKE_CXX_STANDARD="11" \
	-DCMAKE_Fortran_COMPILER=${MPI_HOME}/bin/mpif77 \
	-DCMAKE_Fortran_FLAGS="-fPIC -fallow-argument-mismatch -w" \
	-DBoost_NO_BOOST_CMAKE=TRUE \
	-DHAVE_ACRO=OFF \
	-DHAVE_JEGA=OFF \
	-C${DAK_SRC}/cmake/BuildDakotaCustom.cmake \
	-C${DAK_SRC}/cmake/DakotaDev.cmake \
	${DAK_SRC}

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

cd ${DAK_INSTALL}

# Comment out definition of HAVE_MPI in Teuchos config header file in order to
# avoid conflict with our definition
sed -i -e "s/#define HAVE_MPI/\/* #define HAVE_MPI *\//g" include/Teuchos_config.h
