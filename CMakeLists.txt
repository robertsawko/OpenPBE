cmake_minimum_required (VERSION 2.6)
project (PBESolvers CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#The build is based on OpenFOAM environmental varialbles Need to check 
#if they are loaded and if not TODO search and load an appropriate file?
if(DEFINED ENV{WM_PROJECT})
  message("\nOpenFOAM $ENV{WM_PROJECT_VERSION} found in $ENV{WM_PROJECT_DIR}.\n")
else()
  message("WM_PROJECT is not defined. You must set OpenFOAM specific environment
  variables for CMake to configure the build.")
  # exit early 
  return()
endif()

find_package(Eigen3 REQUIRED)
include_directories (${EIGEN3_INCLUDE_DIR})

#TODO Is there a better way to break lines?
set(CMAKE_CXX_FLAGS "-g -m64 -Dlinux64 -DWM_DP") 
#g++ will inform you of unused typedefs and we want to suppress this as there
#are typedefs like this in OF

if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-local-typedefs")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNoRepository -ftemplate-depth-100")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Xlinker --add-needed -Xlinker --no-as-needed")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_HAS_TR1_TUPLE=0 -wd327,654,819,1125,1476,1505,1572 -xSSE3 -O2 -no-prec-div  -DNoRepository")
endif()

#define a macro to simplify adding OpenFOAM libraries
#TODO  add for each
macro (TARGET_LINK_OPENFOAM_LIBRARIES target)
  foreach(arg ${ARGN})
    if(NOT DEFINED ${arg}_LIB)
      find_library(${arg}_LIB lib${arg}.so $ENV{FOAM_LIBBIN})
    endif()
    target_link_libraries(${target} ${${arg}_LIB})
  endforeach(arg ${ARGN})
endmacro (TARGET_LINK_OPENFOAM_LIBRARIES)

#TODO  add for each
macro (TARGET_LINK_OPENFOAM_USERLIB arg1 arg2)
  if(NOT DEFINED ${arg2}_LIB)
    find_library(${arg2}_LIB lib${arg2}.so $ENV{FOAM_USER_LIBBIN})
  endif()
  #message("Using ${${arg2}_LIB}")
  target_link_libraries(${arg1} ${${arg2}_LIB})
endmacro (TARGET_LINK_OPENFOAM_USERLIB)

enable_testing()

#option (UseOldLibs "Use old libraries" OFF)
#General includes required by all code
#include_directories ("${PROJECT_SOURCE_DIR}/src/multiphaseTurbulence/turbulenceModel")
include_directories ("$ENV{FOAM_SRC}/OpenFOAM/lnInclude")
include_directories ("$ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude")

#Less general but seems to be used by everything too?
include_directories ("$ENV{FOAM_SRC}/transportModels")
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")

#We need to tell CMake to use use C files as CXX (OpenFOAM standard)
set_source_files_properties(filename.c PROPERTIES LANGUAGE CXX )

#gtest for testing
add_subdirectory(ext/gtest/)
#find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
find_package(Boost REQUIRED)

# add_subdirectory(ext/eigen3/)

# add_subdirectory(ext/eigen3/)
#Library
add_subdirectory(src/PBESystem)
add_subdirectory(applications/PBEFoam)


#unit tests
find_package (Threads)
add_subdirectory(tests)

#benchmarks
set(CMAKE_BUILD_TYPE Release)
add_subdirectory(ext/benchmark)
set(GBENCH_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/ext/benchmark/include)
add_subdirectory(benchmarks)
