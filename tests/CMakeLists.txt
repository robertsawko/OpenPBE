set(test_SRCS
    breakup.cpp
    integration.cpp
)

include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories(${GTEST_INCLUDE_DIRS})
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")
add_executable(unit_tests ${test_SRCS})

# target_link_openfoam_userlib(unit_tests PBE)
target_link_libraries(unit_tests PBE)
target_link_openfoam_libraries(unit_tests compressibleTwoPhaseSystem
  compressibleEulerianInterfacialModels compressibleTurbulenceModels
  incompressibleTurbulenceModels phaseCompressibleTurbulenceModels
fvOptions)

target_link_libraries(unit_tests ${CMAKE_BINARY_DIR}/ext/gtest/src/googletest-build/libgtest_main.a)
target_link_libraries(unit_tests ${CMAKE_BINARY_DIR}/ext/gtest/src/googletest-build/libgtest.a)
target_link_libraries(unit_tests ${CMAKE_THREAD_LIBS_INIT})

add_test(binary_breakup_value_test unit_tests)

