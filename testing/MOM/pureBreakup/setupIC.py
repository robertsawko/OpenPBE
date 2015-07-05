from os import path
from numpy import arange
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from scipy.stats import gamma

class_number = 80
vmax = 20.0
dv = vmax / class_number
m0 = 2
prob = gamma(a=7.5, scale=1.0)
v = dv + dv * arange(class_number)

if __name__ == "__main__":
    template_case = SolutionDirectory(
        "MOC-template", archive=None, paraviewLink=False)
    case = template_case.cloneCase(
        "{0}{1}".format("MOC", class_number)
    )

    n0 = ParsedParameterFile(path.join(template_case.name, "0", "n0"))
    for i in range(class_number):
        n0.header["object"] = "n" + str(i)
        n0["internalField"].setUniform(
            m0 * (prob.cdf(v[i]) - prob.cdf(v[i] - dv)) / dv
        )
        n0.writeFileAs(path.join(case.name, "0", "n" + str(i)))

    phase_properties = ParsedParameterFile(
        path.join(case.name, "constant", "phaseProperties")
    )

    phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["numberOfClasses"] = class_number
    phase_properties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["xi1"] = dv
    phase_properties["blending"]["default"]["type"] = "none"
    phase_properties["drag"][1]["swarmCorrection"]["type"] = "none"
    phase_properties.writeFile()

    control_dict = ParsedParameterFile(
        path.join(case.name, "system", "controlDict")
    )
    control_dict["functions"]["probes"]["fields"] = [
        "n{0}".format(n) for n in range(class_number)
    ]
    control_dict.writeFile()
