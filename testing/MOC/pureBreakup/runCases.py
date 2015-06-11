# from os import path
from PyFoam.Execution.BasicRunner import BasicRunner

from parameterizedVariation import nr_classes

if __name__ == "__main__":

    for N in nr_classes:
        block_mesh = BasicRunner(
            argv=["blockMesh", "-case", "testCase{0}".format(N)],
            silent=True
        )
        block_mesh.start()
        run = BasicRunner(
            argv=["twoPhaseEulerFoam", "-case", "testCase{0}".format(N)],
            silent=True
        )
        run.start()
