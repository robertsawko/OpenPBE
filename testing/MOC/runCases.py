# from os import path
from PyFoam.Execution.BasicRunner import BasicRunner

from setup_cases import pbe_grids

if __name__ == "__main__":

    for c in ["pure_coal", "pure_br"]:
        for N in pbe_grids:
            block_mesh = BasicRunner(
                argv=["blockMesh", "-case", "{0}{1}".format(c, N)],
                silent=True
            )
            block_mesh.start()
            run = BasicRunner(
                argv=["twoPhaseEulerFoam", "-case", "{0}{1}".format(c, N)],
                silent=True
            )
            run.start()
