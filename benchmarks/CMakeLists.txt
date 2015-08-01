set(benchmark_SRCS
    moment_inversion.cpp
    integration.cpp
    breakup.cpp
    coalescence.cpp
    main.cpp
)

include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories(${GBENCH_INCLUDE_DIRS})
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")
include_directories ("${PROJECT_SOURCE_DIR}/src/PBESystem")
add_executable(benchmarks ${benchmark_SRCS})

target_link_libraries(benchmarks PBE)
target_link_openfoam_libraries(benchmarks compressibleTwoPhaseSystem
  compressibleEulerianInterfacialModels compressibleTurbulenceModels
  incompressibleTurbulenceModels phaseCompressibleTurbulenceModels
fvOptions)

target_link_libraries(benchmarks ${CMAKE_BINARY_DIR}/ext/benchmark/src/libbenchmark.a)
target_link_libraries(benchmarks ${CMAKE_THREAD_LIBS_INIT})

