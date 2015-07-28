from os import path
from numpy import exp, arange
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from generate_initial_conditions import calculate_MOC_initial_condition


class test_case:
    def __init__(
        self, nr_classes
    ):
        self.name = "re6400"
        self.nr_classes = nr_classes
        self.dv, self.Ninit = calculate_MOC_initial_condition(nr_classes)


def case_setup(ci):
    template_case = SolutionDirectory(
        "template", archive=None, paraviewLink=False)
    case = template_case.cloneCase(
        "{0}NC{1}".format(ci.name, ci.nr_classes)
    )

    phase_properties = ParsedParameterFile(
        path.join(template_case.name, "constant", "phaseProperties"))
    phase_properties["oil"]["PBEDiameterCoeffs"]["MOCCoeffs"]["numberOfClasses"] = ci.nr_classes
    phase_properties["oil"]["PBEDiameterCoeffs"]["MOCCoeffs"]["xi1"] = ci.dv

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

    controlDict = ParsedParameterFile(
        path.join(case.name, "system", "controlDict")
    )
    controlDict["functions"]["probes"]["fields"] = [
        "n{0}".format(n) for n in range(ci.nr_classes)]
    controlDict.writeFile()

pbe_grids = [10, 20, 40, 80]
galinat_cases = [test_case(g) for g in pbe_grids]

if __name__ == "__main__":

    for c in galinat_cases:
        case_setup(c)
