# MacBook
set(CMAKE_C_COMPILER   "cc")
set(CMAKE_CXX_COMPILER "c++")
set (CXX_COMPILER_WRAPPER mpicxx)
set (C_COMPILER_WRAPPER mpicc)

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(FFTW_INCLUDE_DIR   "/opt/local/include")
set(FFTW_LIB           "/opt/local/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/opt/local/include")
set(NETCDF_LIB_C       "/opt/local/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/opt/local/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/opt/local/lib/libhdf5.a")
set(HDF5_LIB_2         "/opt/local/lib/libhdf5_hl.a")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
