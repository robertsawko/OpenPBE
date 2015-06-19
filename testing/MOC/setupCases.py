from os import path
from numpy import exp, arange
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

nr_classes = [10, 20, 40, 80]


class test_case:
    def __init__(
        self, case_name, nr_classes, vmax,
        N0, phase_properties_name,
        end_time, delta_t
    ):
        self.dv = vmax / nC
        template_case = SolutionDirectory(
            "template", archive=None, paraviewLink=False)
        self.case = template_case.cloneCase("{0}{1}".format(case_name, nr_classes))

        phase_properties = ParsedParameterFile(
            path.join("./diffs", phase_properties_name))
        phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["numberOfClasses"] = nr_classes
        phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["xi1"] =\
            self.dv

        # manually fix bad pyfoam parsing
        phase_properties["blending"]["default"]["type"] = "none"
        phase_properties["drag"][1]["swarmCorrection"]["type"] = "none"
        phase_properties.writeFileAs(path.join(
            self.case.name, "constant", "phaseProperties"
        ))

        v = self.dv + self.dv * arange(nr_classes)
        n0 = ParsedParameterFile(path.join(template_case.name, "0", "n0"))
        for i in range(nr_classes):
            n0.header["object"] = "n" + str(i)
            n0["internalField"].setUniform(N0(v[i]))
            n0.writeFileAs(path.join(self.case.name, "0", "n" + str(i)))

        controlDict = ParsedParameterFile(
            path.join(self.case.name, "system", "controlDict")
        )
        controlDict["functions"]["probes"]["fields"] = [
            "n{0}".format(n) for n in range(nr_classes)]
        controlDict["endTime"] = end_time
        controlDict["deltaT"] = delta_t
        controlDict.writeFile()


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
        self.vmax = 1e1
        # C = 0.1
        N0 = 2
        v0 = 0.5
        test_case.__init__(
            self, "pure_coal", nr_classes, self.vmax,
            lambda v: (N0 / v0) * (v / v0) * exp(-v / v0) * self.dv,
            "phaseProperties.coalescence",
            1, 0.001
        )

if __name__ == "__main__":

    for nC in nr_classes:
        coalescence_case(nC)
        breakup_case(nC)
