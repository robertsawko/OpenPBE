include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories(${GTEST_INCLUDE_DIRS})
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")
include_directories ("${PROJECT_SOURCE_DIR}/src/PBESystem")

function (ADD_UNIT_TEST testname executable sourcefile)

  set(sources ${sourcefile})
  add_executable(${executable} ${sources})

  target_link_libraries(${executable} PBE)
  target_link_openfoam_libraries(${executable} compressibleTwoPhaseSystem
    compressibleEulerianInterfacialModels compressibleTurbulenceModels
    incompressibleTurbulenceModels phaseCompressibleTurbulenceModels
    fvOptions)

  target_link_libraries(${executable} ${CMAKE_BINARY_DIR}/ext/gtest/src/googletest-build/libgtest_main.a)
  target_link_libraries(${executable} ${CMAKE_BINARY_DIR}/ext/gtest/src/googletest-build/libgtest.a)
  target_link_libraries(${executable} ${CMAKE_THREAD_LIBS_INIT})

  add_test(${testname} ${executable})
endfunction (ADD_UNIT_TEST)

add_unit_test(breakup_unittest breakup_test breakup.cpp integration.cpp)
add_unit_test(numerical_integration_unittest integration_test integration.cpp)
add_unit_test(moment_inversion_unittestest moment_inversion_test moment_inversion.cpp)


