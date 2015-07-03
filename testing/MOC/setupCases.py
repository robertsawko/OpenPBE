from os import path
from numpy import exp, arange
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

pbe_grids = [10, 20, 40, 80]


class test_case:
    def __init__(
        self, case_name, nr_classes, vmax,
        N0, phase_properties_name,
        end_time, delta_t
    ):
        self.name = case_name
        self.nr_classes = nr_classes
        self.vmax = vmax
        self.N0 = N0
        self.phase_properties_name = phase_properties_name
        self.end_time = end_time
        self.delta_t = delta_t


class breakup_case(test_case):

    def __init__(self, nr_classes):
        self.l = 1.0

        def Nbr0(v):
            if(v < self.l):
                return 0
            else:
                return 1

        test_case.__init__(
            self, "pure_br", nr_classes, self.l,
            Nbr0,
            "phaseProperties.breakup",
            10, 0.005

        )


class coalescence_case(test_case):

    def __init__(self, nr_classes):
        # C = 0.1
        vmax = 1e1
        dv = float(vmax) / nr_classes
        N0 = 2
        v0 = 0.5
        test_case.__init__(
            self, "pure_coal", nr_classes, vmax,
            lambda v: (N0 / v0) * (v / v0) * exp(-v / v0) * dv,
            "phaseProperties.coalescence",
            1, 0.001
        )


def case_setup(ci):
    dv = ci.vmax / ci.nr_classes
    template_case = SolutionDirectory(
        "template", archive=None, paraviewLink=False)
    case = template_case.cloneCase(
        "{0}{1}".format(ci.name, ci.nr_classes)
    )

    phase_properties = ParsedParameterFile(
        path.join("./diffs", ci.phase_properties_name))
    phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["numberOfClasses"] = ci.nr_classes
    phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["xi1"] = dv

    # manually fix bad pyfoam parsing
    phase_properties["blending"]["default"]["type"] = "none"
    phase_properties["drag"][1]["swarmCorrection"]["type"] = "none"
    phase_properties.writeFileAs(path.join(
        case.name, "constant", "phaseProperties"
    ))

    v = dv + dv * arange(ci.nr_classes)
    n0 = ParsedParameterFile(path.join(template_case.name, "0", "n0"))
    for i in range(ci.nr_classes):
        n0.header["object"] = "n" + str(i)
        n0["internalField"].setUniform(ci.N0(v[i]))
        n0.writeFileAs(path.join(case.name, "0", "n" + str(i)))

    controlDict = ParsedParameterFile(
        path.join(case.name, "system", "controlDict")
    )
    controlDict["functions"]["probes"]["fields"] = [
        "n{0}".format(n) for n in range(ci.nr_classes)]
    controlDict["endTime"] = ci.end_time
    controlDict["deltaT"] = ci.delta_t
    controlDict.writeFile()


if __name__ == "__main__":

    for nC in pbe_grids:
        case_setup(coalescence_case(nC))
        case_setup(breakup_case(nC))
