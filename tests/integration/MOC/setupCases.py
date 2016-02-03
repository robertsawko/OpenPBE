from os import path
from numpy import exp, arange
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


class test_case:
    def __init__(
        self, case_name, nr_classes, vmax,
        Ninit, phase_properties_name,
        end_time, delta_t
    ):
        self.name = case_name
        self.nr_classes = nr_classes
        self.vmax = vmax
        self.dv = float(vmax) / nr_classes
        self.Ninit = Ninit
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
            self, "pure_br", nr_classes, vmax=self.l,
            Ninit=Nbr0, phase_properties_name="phaseProperties.breakup",
            end_time=10, delta_t=0.005
        )


class coalescence_case(test_case):

    def __init__(self, nr_classes, N0=2, v0=0.5, C=0.1):
        vmax = 1e1
        dv = float(vmax) / nr_classes
        self.N0 = N0
        self.v0 = v0
        self.C = C
        test_case.__init__(
            self, "pure_coal", nr_classes, vmax,
            Ninit=lambda v: (N0 / v0) * (v / v0) * exp(-v / v0) * dv,
            phase_properties_name="phaseProperties.coalescence",
            end_time=1, delta_t=0.001
        )


def case_setup(ci):
    template_case = SolutionDirectory(
        "template", archive=None, paraviewLink=False)
    case = template_case.cloneCase(
        "{0}{1}".format(ci.name, ci.nr_classes)
    )

    phase_properties = ParsedParameterFile(
        path.join("./diffs", ci.phase_properties_name))
    phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["numberOfClasses"] = ci.nr_classes
    phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["xi1"] = ci.dv

    # manually fix bad pyfoam parsing
    phase_properties["blending"]["default"]["type"] = "none"
    phase_properties["drag"][1]["swarmCorrection"]["type"] = "none"
    phase_properties.writeFileAs(path.join(
        case.name, "constant", "phaseProperties"
    ))

    v = ci.dv + ci.dv * arange(ci.nr_classes)
    n0 = ParsedParameterFile(path.join(template_case.name, "0", "n0"))
    for i in range(ci.nr_classes):
        n0.header["object"] = "n" + str(i)
        n0["internalField"].setUniform(ci.Ninit(v[i]))
        n0.writeFileAs(path.join(case.name, "0", "n" + str(i)))

    controlDict = ParsedParameterFile(
        path.join(case.name, "system", "controlDict")
    )
    controlDict["functions"]["probes"]["fields"] = [
        "n{0}".format(n) for n in range(ci.nr_classes)]
    controlDict["endTime"] = ci.end_time
    controlDict["deltaT"] = ci.delta_t
    controlDict.writeFile()

pbe_grids = [10, 20, 40, 80]
coalescence_cases = [coalescence_case(nC) for nC in pbe_grids]
breakup_cases = [breakup_case(nC) for nC in pbe_grids]

if __name__ == "__main__":

    for bc in breakup_cases:
        case_setup(bc)

    for cc in coalescence_cases:
        case_setup(cc)
