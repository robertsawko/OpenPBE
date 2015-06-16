# find_package(Eigen3 REQUIRED)

include_directories (EIGEN3_INCLUDE_DIRS)
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")

set(lib_header_files
  coalescenceKernels/coalescenceKernel.H
  coalescenceKernels/noCoalescence.H
  breakupKernels/binaryBreakup.H
  breakupKernels/breakupKernel.H
  breakupKernels/noBreakup.H
  daughterParticleDistributions/daughterParticleDistribution.H
  daughterParticleDistributions/uniformBinaryBreakup.H
  PBEDiameter.H
  PBEMethod.H
  dummy.H
  MOC.H
  MOM.H
  PBESystems-internal.C
)

set(lib_source_files
  coalescenceKernels/coalescenceKernel.C
  coalescenceKernels/noCoalescence.C
  breakupKernels/binaryBreakup.C
  breakupKernels/breakupKernel.C
  breakupKernels/noBreakup.C
  daughterParticleDistributions/daughterParticleDistribution.C
  daughterParticleDistributions/uniformBinaryBreakup.C
  PBEDiameter.C
  PBEMethod.C
  dummy.C
  MOC.C
  MOM.C
)

add_library(PBE SHARED ${lib_header_files} ${lib_source_files})

set_target_properties (PBE PROPERTIES LIBRARY_OUTPUT_DIRECTORY
  "$ENV{FOAM_USER_LIBBIN}")
