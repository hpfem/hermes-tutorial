project(hermes-tutorial)
	cmake_minimum_required(VERSION 2.8)
  
###########################################################################
#############//////---------------------------------\\\\\\#################
#############!!!!!! NEVER CHANGE THIS FILE DIRECTLY !!!!!!#################
################!!!!!! USE the "CMake.vars" file !!!!!!####################
#############\\\\\\---------------------------------//////#################
###########################################################################
#   __      __      __      __      __      __      __      __      __    #
# /    \  /    \  /    \  /    \  /    \  /    \  /    \  /    \  /    \  #
#| STOP || STOP || STOP || STOP || STOP || STOP || STOP || STOP || STOP | #
# \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#  /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \   #
# ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ #
###########################################################################
#############//////---------------------------------\\\\\\#################
#############!!!!!! NEVER CHANGE THIS FILE DIRECTLY !!!!!!#################
################!!!!!! USE the "CMake.vars" file !!!!!!####################
#############\\\\\\---------------------------------//////#################
###########################################################################

  set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
  include(CommonTargetProperties)

  # For Win64
  if(${CMAKE_CL_64})
    set(WIN64 YES)
  else(${CMAKE_CL_64})
    set(WIN64 NO)
  endif(${CMAKE_CL_64})

  # Where to look for the static libraries.
  set(HERMES_DIRECTORY /usr/local/lib)
  set(HERMES_INCLUDE_PATH /usr/local/include)
  set(HERMES_SOURCE_DIRECTORY ~/hermes)

  # Allow to override the default values in CMake.vars:
  include(CMake.vars OPTIONAL)
  include(".CMake.vars.default.Linux" OPTIONAL)
  include(${HERMES_CMAKE_VARS_FILE_LOCATION} OPTIONAL)

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		# Disable all warnings and turn on only important ones:
		set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} -std=c++11")
		set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")

		set(RELEASE_FLAGS "-DNDEBUG -O3")
		set(DEBUG_FLAGS   "-g")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  # This overrides CXX flags for MSVC
  if(MSVC)
    include(MSVC_FLAGS)
  endif(MSVC)

  include_directories("${PTHREAD_ROOT}/include")

  if(${WITH_UMFPACK})
    find_package(UMFPACK REQUIRED)
    include_directories(${UMFPACK_INCLUDE_DIRS})
  endif()

  find_package(PTHREAD REQUIRED)

  if(${WITH_TC_MALLOC})
    find_package(TCMALLOC REQUIRED)
    include_directories(${TCMALLOC_INCLUDE_DIR})
  endif()

  if(${WITH_PARALUTION})
    find_package(PARALUTION REQUIRED)
    include_directories(${PARALUTION_INCLUDE_DIR})
  endif()

  if(${WITH_MATIO})
    find_package(MATIO REQUIRED)
    include_directories(${MATIO_INCLUDE_DIR})
  endif()

  if(MSVC)
    if(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLU)
        find_package(CLAPACK REQUIRED)
        include_directories(${CLAPACK_INCLUDE_DIRS})
    endif(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLU)
    if(WITH_MUMPS)
      find_package(WINBLAS REQUIRED)
    endif(WITH_MUMPS)
  else(MSVC)
    if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
      if (NOT LAPACK_FOUND)
          enable_language(Fortran)
          find_package(LAPACK REQUIRED)
          set(LAPACK_LIBRARIES  ${LAPACK_LIBRARIES}
                  CACHE STRING  "Path to LAPACK/BLAS libraries.")
          set(LAPACK_FOUND      YES
                  CACHE STRING  "Have LAPACK/BLAS libraries been found?")
      endif (NOT LAPACK_FOUND)
      add_definitions(-DWITH_BLAS)  
    endif(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
  endif(MSVC)

  if(WITH_MUMPS)
    find_package(MUMPS REQUIRED)
    if(WITH_MPI)
      find_package(SCALAPACK REQUIRED) # parallel version of MUMPS needs ScaLapack+BLACS
    endif(WITH_MPI)
    include_directories(${MUMPS_INCLUDE_DIR})
    set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES})
    LIST(APPEND MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES})
  endif(WITH_MUMPS)

  if(WITH_TRILINOS)
    find_package(TRILINOS REQUIRED)
    include_directories(${TRILINOS_INCLUDE_DIR})
  endif(WITH_TRILINOS)
      
  if(${WITH_BSON})
    find_package(BSON REQUIRED)
    include_directories(${BSON_INCLUDE_DIR})
  endif()
        
  if(WITH_SUPERLU)
    find_package(SUPERLU REQUIRED)
    include_directories(${SUPERLU_INCLUDE_DIR})
    if(SUPERLU_MT)
      add_definitions(-DSLU_MT)
    endif(SUPERLU_MT)
  endif(WITH_SUPERLU)
  set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
  set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

  find_package(HERMES_COMMON REQUIRED)

  include_directories(${HERMES_COMMON_INCLUDE_PATH})
  include_directories(${HERMES2D_INCLUDE_PATH})

  # Linear problems.
  add_subdirectory(A-linear)

  # Nonlinear problems.
  add_subdirectory(B-nonlinear)

  # Time-dependent problems.
  add_subdirectory(C-transient)

  # Adaptivity.
  add_subdirectory(D-adaptivity)

  # FVM and DG.
  add_subdirectory(E-fvm-and-dg)

  # Using Trilinos.
  if(WITH_TRILINOS)
    add_subdirectory(F-trilinos)
  endif(WITH_TRILINOS)

  # Miscellaneous techniques.
  add_subdirectory(G-miscellaneous)

  #--------
  # Report.
  #~~~~~~~~

  message("\nHermes Tutorial Configuration results")
  message("-------------------------------------")
  if(WIN32)
    message("Build for Win64: ${WIN64}")
  endif(WIN32)
  
  message("\n--------------Features---------------")
  message("Build with OpenMP: ${WITH_OPENMP}")
  message("Build with TCMalloc: ${WITH_TC_MALLOC}")
  message("Build with BSON: ${WITH_BSON}")
  message("Build with MATIO: ${WITH_MATIO}")
  if(${WITH_MATIO})
    message(" MATIO with HDF5: ${MATIO_WITH_HDF5}")
  endif()
  if(${WITH_MPI})
    message("Build with MPI: ${WITH_MPI}")
  endif()
  message("Build with EXODUSII: ${WITH_EXODUSII}")
  message("-------------------------------------")
  
  message("\n---------------Solvers---------------")
  message("Build with UMFPACK: ${WITH_UMFPACK}")
  message("Build with PARALUTION: ${WITH_PARALUTION}")
  message("Build with PETSC: ${WITH_PETSC}")
  message("Build with MUMPS: ${WITH_MUMPS}")
  message("Build with SUPERLU${MT}: ${WITH_SUPERLU}")
  message("Build with TRILINOS: ${WITH_TRILINOS}")
  message("-------------------------------------")
  
  if(WITH_H2D)
    message("\n--------------Hermes2D---------------")
    message(" Hermes2D with OpenGL: ${H2D_WITH_GLUT}")
  endif(WITH_H2D)
  message("-------------------------------------")
  message("\n")