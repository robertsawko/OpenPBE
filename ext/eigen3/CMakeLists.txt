cmake_minimum_required(VERSION 2.8.8)
project(eigen_checkout C CXX)
include(ExternalProject)

ExternalProject_Add(eigen
    HG_REPOSITORY https://bitbucket.org/eigen/eigen/
    HG_TAG 5a34505
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
    # Disable configure, build, install (only checkout)
    CONFIGURE_COMMAND ""
    CMAKE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Specify include dir
ExternalProject_Get_Property(eigen source_dir)
set(EIGEN3_INCLUDE_DIRS ${source_dir}/include PARENT_SCOPE)
