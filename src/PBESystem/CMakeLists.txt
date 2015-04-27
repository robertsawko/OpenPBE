# file (GLOB test_sources
#   MOMTests.C
#   PBESystems-internal.C
# )

# add_executable(runMOMTests ${test_sources})
# Link test executable against gtest
# add_dependencies(runMOMTests googletest)
# target_link_libraries(runMOMTests ${GTEST_LIBS_DIR}/libgtest.a)
# target_link_openfoam_libraries(runMOMTests finiteVolume
#     fvOptions OpenFOAM) 


file(GLOB lib_source_files
  coalescenceKernels/coalescenceKernel.C
  coalescenceKernels/noCoalescence.C
  breakupKernels/binaryBreakup.C
  breakupKernels/breakupKernel.C
  breakupKernels/noBreakup.C
  PBEDiameter.C
  PBEMethod.C
  dummy.C
  MOC.C
)

# file(GLOB test_files
#   *Tests.C
# )
# list(REMOVE_ITEM lib_source_files 
#   ${test_files})


add_library(PBE SHARED ${lib_source_files})
include_directories("${PROJECT_SOURCE_DIR}/src/eigen")

#include_directories ("$ENV{FOAM_APP}/solvers/multiphase/multiphaseEulerFoam/phaseModel/lnInclude")
#include_directories ("PBEPhaseModel/phaseModel")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
#include_directories ("${PROJECT_SOURCE_DIR}/src/multiphaseTurbulence/")
include_directories ("./")

set_target_properties (PBE PROPERTIES LIBRARY_OUTPUT_DIRECTORY
  "$ENV{FOAM_USER_LIBBIN}")


# add_test(MOMTests runMOMTests)
