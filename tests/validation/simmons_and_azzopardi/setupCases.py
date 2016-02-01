from os import path
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from numpy import zeros, pi


class test_case:
    def __init__(
        self, nr_classes
    ):
        self.name = "c01"
        self.nr_classes = nr_classes
        dmax = 7.0e-4  # Read from figure9 (could calculate from Hinze!)
        vmax = pi / 6 * dmax**3
        self.dv, self.Ninit = vmax / nr_classes, zeros(nr_classes)
        self.Ninit[nr_classes // 2] = 1


def case_setup(ci):
    template_case = SolutionDirectory(
        "c01", archive=None, paraviewLink=False)
    case = template_case.cloneCase(
        "{0}NC{1}".format(ci.name, ci.nr_classes)
    )

    phase_properties = ParsedParameterFile(
        path.join(template_case.name, "constant", "phaseProperties"))
    a, b, c = "acps", "PBEDiameterCoeffs", "MOCCoeffs"
    phase_properties[a][b][c]["numberOfClasses"] = ci.nr_classes
    phase_properties[a][b][c]["xi1"] = ci.dv

    # manually fix bad pyfoam parsing
    phase_properties["blending"]["default"]["type"] = "none"
    phase_properties["drag"][1]["swarmCorrection"]["type"] = "none"
    phase_properties.writeFileAs(path.join(
        case.name, "constant", "phaseProperties"
    ))

    n0 = ParsedParameterFile(path.join(template_case.name, "0", "n0"))
    for i in range(ci.nr_classes):
        n0.header["object"] = "n" + str(i)
        n0["internalField"].setUniform(ci.Ninit[i])
        n0.writeFileAs(path.join(case.name, "0", "n" + str(i)))

    # controlDict = ParsedParameterFile(
        # path.join(case.name, "system", "controlDict")
    # )
    # controlDict["functions"]["probes"]["fields"] = [
    #     "n{0}".format(n) for n in range(ci.nr_classes)]
    # controlDict.writeFile()

pbe_grids = [50]
cases = [test_case(g) for g in pbe_grids]

if __name__ == "__main__":

    for c in cases:
        case_setup(c)
