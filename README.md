#OpenPBE


## Objectives

 * Develop the PBE port for the new solver.
 * Port Marcin's solution method from of2.2.x.
 * Implement method of classes (MOC)
 * Implement quadrature method of moments (QMOM)

## Build

```bash
# Clone the project

# Load OpenFOAM 2.3.1 environment

cd <project_dir>

# Donload benchmark
git submodule init
git submodule update 

# Build the project
mkdir build
cmake ../
make #alternatively make -j 4 for parallell build
```

## Usage

The OF2.3 branch uses `twoPhaseEulerFoam` from the main distribution. The only
necessary change is to include

```
libs ("libPBE.so")
``` 

in your `controlDict`. See `testing/twoPhaseEulerFoam` for the full setup.
