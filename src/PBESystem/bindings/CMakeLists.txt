find_package(PythonLibs 2)
find_package(Boost COMPONENTS python REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_library(kernels SHARED kernels.cpp)
set_target_properties(kernels PROPERTIES PREFIX "")
target_link_libraries (kernels
    boost_python
    ${PYTHON_LIBRARIES}
    ${Boost_LIBRARIES}
    PBE
)
target_link_openfoam_libraries(kernels compressibleTwoPhaseSystem
  compressibleEulerianInterfacialModels compressibleTurbulenceModels
  incompressibleTurbulenceModels phaseCompressibleTurbulenceModels
fvOptions)
