# SuperMUC
set(CMAKE_C_COMPILER   "icc")
set(CMAKE_CXX_COMPILER "icc")
set (CXX_COMPILER_WRAPPER mpicxx)
set (C_COMPILER_WRAPPER mpicc)

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(FFTW_INCLUDE_DIR   "/lrz/sys/libraries/fftw/3.3.2/sse/include")
set(FFTW_LIB           "/lrz/sys/libraries/fftw/3.3.2/sse/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/lrz/sys/libraries/netcdf/4.2.1.1/include")
set(NETCDF_LIB_C       "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5.a")
set(HDF5_LIB_2         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5_hl.a")
set(SZIP_LIB           "/lrz/sys/libraries/hdf5/szip_2.1_u1/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
