add_executable(PBEFoam PBEFoam.C)

#Includes specfic to the project

#General CFD libraries
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseCompressible/lnInclude")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/phaseCompressibleTurbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/interfacialModels/lnInclude")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")

include_directories ("$ENV{FOAM_SRC}/fvOptions/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")
include_directories ("$ENV{FOAM_SRC}/sampling/lnInclude")

#Solver specific libraries
include_directories ("${PROJECT_SOURCE_DIR}/src/PBESystem")

target_link_openfoam_libraries(PBEFoam  
    fluidThermophysicalModels specie turbulenceModels
    compressibleTurbulenceModels phaseCompressibleTurbulenceModels
    incompressibleTransportModels compressibleTwoPhaseSystem
    compressibleEulerianInterfacialModels finiteVolume fvOptions meshTools
    sampling OpenFOAM)

target_link_libraries(PBEFoam PBE)

set_target_properties(PBEFoam PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$ENV{FOAM_USER_APPBIN}")

