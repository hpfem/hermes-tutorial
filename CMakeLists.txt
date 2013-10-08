project(hermes-tutorialMeshSharedPtr mesh
cmake_minimum_required(VERSION 2.8MeshSharedPtr mesh

set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmakeMeshSharedPtr mesh
include(CommonTargetPropertiesMeshSharedPtr mesh

# For Win64
if(${CMAKE_CL_64}MeshSharedPtr mesh
  set(WIN64 YESMeshSharedPtr mesh
else(${CMAKE_CL_64}MeshSharedPtr mesh
  set(WIN64 NOMeshSharedPtr mesh
endif(${CMAKE_CL_64}MeshSharedPtr mesh

# This has to be the same as in the library's CMake.
# Default "yes"
set(H2D_WITH_GLUT YESMeshSharedPtr mesh
set(WITH_TC_MALLOC YESMeshSharedPtr mesh
set(WITH_PARALUTION YESMeshSharedPtr mesh
set(WITH_UMFPACK YESMeshSharedPtr mesh
set(WITH_OPENMP YESMeshSharedPtr mesh
# Default "no"
set(WITH_EXODUSII NOMeshSharedPtr mesh
set(WITH_MUMPS NOMeshSharedPtr mesh
set(WITH_MATIO NOMeshSharedPtr mesh
set(WITH_BSON NOMeshSharedPtr mesh
set(WITH_SUPERLU NOMeshSharedPtr mesh

# Where to look for the static libraries.
set(HERMES_DIRECTORY /usr/local/libMeshSharedPtr mesh
set(HERMES_INCLUDE_PATH /usr/local/includeMeshSharedPtr mesh
set(DEP_INCLUDE_PATHS /usr/local/includeMeshSharedPtr mesh

# Allow to override the default values in CMake.vars:
include(CMake.vars OPTIONALMeshSharedPtr mesh

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"MeshSharedPtr mesh
		# Disable all warnings and turn on only important ones:
		set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} -std=c++11"MeshSharedPtr mesh
		set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}"MeshSharedPtr mesh

		set(RELEASE_FLAGS "-DNDEBUG -O3"MeshSharedPtr mesh
		set(DEBUG_FLAGS   "-g"MeshSharedPtr mesh
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"MeshSharedPtr mesh

# Enabling multiprocessor build on MSVC
if(MSVCMeshSharedPtr mesh
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP"MeshSharedPtr mesh
  if(WITH_OPENMPMeshSharedPtr mesh
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp"MeshSharedPtr mesh
  endif(WITH_OPENMPMeshSharedPtr mesh
endif(MSVCMeshSharedPtr mesh

# This overrides CXX flags for MSVC
if(MSVCMeshSharedPtr mesh
  if(WIN64MeshSharedPtr mesh
    set(MSVC_DEFINES "/DWIN64 /D_WINDOWS /Dfinite=_finite /wd4275 /wd4251"MeshSharedPtr mesh
  else(WIN64MeshSharedPtr mesh
    set(MSVC_DEFINES "/DWIN32 /D_WINDOWS /Dfinite=_finite /wd4275 /wd4251"MeshSharedPtr mesh
  endif(WIN64MeshSharedPtr mesh
  set(CMAKE_CXX_FLAGS_DEBUG          "/D_DEBUG /Od /Ob2 /MDd /Zi ${MSVC_DEFINES}"MeshSharedPtr mesh
  set(HERMES_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}"MeshSharedPtr mesh
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/O2 /Ob2 /MD ${MSVC_DEFINES}"MeshSharedPtr mesh
  set(CMAKE_CXX_FLAGS_RELEASE        "/DNDEBUG /O2 /Ob2 /MD ${MSVC_DEFINES}"MeshSharedPtr mesh
  set(HERMES_RELEASE_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}"MeshSharedPtr mesh
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/DNDEBUG /O2 /Ob2 /MD /Zi ${MSVC_DEFINES}"MeshSharedPtr mesh
  
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO"MeshSharedPtr mesh
  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /SAFESEH:NO"MeshSharedPtr mesh
  set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /SAFESEH:NO"MeshSharedPtr mesh
endif(MSVCMeshSharedPtr mesh

include_directories("${PTHREAD_ROOT}/include"MeshSharedPtr mesh

if(${WITH_UMFPACK}MeshSharedPtr mesh
  find_package(UMFPACK REQUIREDMeshSharedPtr mesh
  include_directories(${UMFPACK_INCLUDE_DIRS}MeshSharedPtr mesh
endif(MeshSharedPtr mesh

find_package(PTHREAD REQUIREDMeshSharedPtr mesh

if(${WITH_TC_MALLOC}MeshSharedPtr mesh
  find_package(TCMALLOC REQUIREDMeshSharedPtr mesh
  include_directories(${TCMALLOC_INCLUDE_DIR}MeshSharedPtr mesh
endif(MeshSharedPtr mesh

if(${WITH_PARALUTION}MeshSharedPtr mesh
  find_package(PARALUTION REQUIREDMeshSharedPtr mesh
  include_directories(${PARALUTION_INCLUDE_DIR}MeshSharedPtr mesh
endif(MeshSharedPtr mesh

if(${WITH_MATIO}MeshSharedPtr mesh
  find_package(MATIO REQUIREDMeshSharedPtr mesh
  include_directories(${MATIO_INCLUDE_DIR}MeshSharedPtr mesh
endif(MeshSharedPtr mesh

if(MSVCMeshSharedPtr mesh
  if(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLUMeshSharedPtr mesh
      find_package(CLAPACK REQUIREDMeshSharedPtr mesh
      include_directories(${CLAPACK_INCLUDE_DIRS}MeshSharedPtr mesh
  endif(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLUMeshSharedPtr mesh
  if(WITH_MUMPSMeshSharedPtr mesh
    find_package(WINBLAS REQUIREDMeshSharedPtr mesh
  endif(WITH_MUMPSMeshSharedPtr mesh
else(MSVCMeshSharedPtr mesh
  if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLUMeshSharedPtr mesh
    if (NOT LAPACK_FOUNDMeshSharedPtr mesh
        enable_language(FortranMeshSharedPtr mesh
        find_package(LAPACK REQUIREDMeshSharedPtr mesh
        set(LAPACK_LIBRARIES  ${LAPACK_LIBRARIES}
                CACHE STRING  "Path to LAPACK/BLAS libraries."MeshSharedPtr mesh
        set(LAPACK_FOUND      YES
                CACHE STRING  "Have LAPACK/BLAS libraries been found?"MeshSharedPtr mesh
    endif (NOT LAPACK_FOUNDMeshSharedPtr mesh
    add_definitions(-DWITH_BLASMeshSharedPtr mesh  
  endif(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLUMeshSharedPtr mesh
endif(MSVCMeshSharedPtr mesh

if(WITH_MUMPSMeshSharedPtr mesh
  find_package(MUMPS REQUIREDMeshSharedPtr mesh
  if(WITH_MPIMeshSharedPtr mesh
    find_package(SCALAPACK REQUIREDMeshSharedPtr mesh # parallel version of MUMPS needs ScaLapack+BLACS
  endif(WITH_MPIMeshSharedPtr mesh
  include_directories(${MUMPS_INCLUDE_DIR}MeshSharedPtr mesh
  set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES}MeshSharedPtr mesh
  LIST(APPEND MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES}MeshSharedPtr mesh
endif(WITH_MUMPSMeshSharedPtr mesh

if(${WITH_BSON}MeshSharedPtr mesh
  find_package(BSON REQUIREDMeshSharedPtr mesh
  include_directories(${BSON_INCLUDE_DIR}MeshSharedPtr mesh
endif(MeshSharedPtr mesh
      
if(WITH_SUPERLUMeshSharedPtr mesh
  find_package(SUPERLU REQUIREDMeshSharedPtr mesh
  include_directories(${SUPERLU_INCLUDE_DIR}MeshSharedPtr mesh
  if(SUPERLU_MTMeshSharedPtr mesh
    add_definitions(-DSLU_MTMeshSharedPtr mesh
  endif(SUPERLU_MTMeshSharedPtr mesh
endif(WITH_SUPERLUMeshSharedPtr mesh
set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2dMeshSharedPtr mesh
set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_commonMeshSharedPtr mesh

find_package(HERMES_COMMON REQUIREDMeshSharedPtr mesh

include_directories(${HERMES_COMMON_INCLUDE_PATH}MeshSharedPtr mesh
include_directories(${HERMES2D_INCLUDE_PATH}MeshSharedPtr mesh
include_directories(${DEP_INCLUDE_PATHS}MeshSharedPtr mesh

# Linear problems.
add_subdirectory(A-linearMeshSharedPtr mesh

# Nonlinear problems.
add_subdirectory(B-nonlinearMeshSharedPtr mesh

# Time-dependent problems.
add_subdirectory(C-transientMeshSharedPtr mesh

# Adaptivity.
add_subdirectory(D-adaptivityMeshSharedPtr mesh

# FVM and DG.
add_subdirectory(E-fvm-and-dgMeshSharedPtr mesh

# Using Trilinos.
if(WITH_TRILINOSMeshSharedPtr mesh
	add_subdirectory(F-trilinosMeshSharedPtr mesh
endif(WITH_TRILINOSMeshSharedPtr mesh

# Miscellaneous techniques.
add_subdirectory(G-miscellaneousMeshSharedPtr mesh

#------------------------------------------------------------------------------
# Report.
#~~~~~~~~
message("\nHermes-tutorial Configuration results"MeshSharedPtr mesh
message("---------------------"MeshSharedPtr mesh
if(WIN32MeshSharedPtr mesh
  message("Build for Win64: ${WIN64}"MeshSharedPtr mesh
endif(WIN32MeshSharedPtr mesh

message("\n-------Features-------"MeshSharedPtr mesh
message("Build with OpenMP: ${WITH_OPENMP}"MeshSharedPtr mesh
message("Build with TCMalloc: ${WITH_TC_MALLOC}"MeshSharedPtr mesh
message("Build with BSON: ${WITH_BSON}"MeshSharedPtr mesh
message("Build with MATIO: ${WITH_MATIO}"MeshSharedPtr mesh
if(${WITH_MATIO}MeshSharedPtr mesh
  message(" MATIO with HDF5: ${MATIO_WITH_HDF5}"MeshSharedPtr mesh
endif(MeshSharedPtr mesh
if(${WITH_MPI}MeshSharedPtr mesh
  message("Build with MPI: ${WITH_MPI}"MeshSharedPtr mesh
endif(MeshSharedPtr mesh
message("Build with EXODUSII: ${WITH_EXODUSII}"MeshSharedPtr mesh

message("\n-------Solvers-------"MeshSharedPtr mesh
message("Build with UMFPACK: ${WITH_UMFPACK}"MeshSharedPtr mesh
message("Build with PARALUTION: ${WITH_PARALUTION}"MeshSharedPtr mesh
message("Build with PETSC: !TO_DO!"MeshSharedPtr mesh
message("Build with MUMPS: ${WITH_MUMPS}"MeshSharedPtr mesh
message("Build with SUPERLU${MT}: ${WITH_SUPERLU}"MeshSharedPtr mesh
message("Build with TRILINOS: !TO_DO!"MeshSharedPtr mesh

message("---------------------"MeshSharedPtr mesh
message("\n"MeshSharedPtr mesh
