from os import path
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from numpy import pi, log, sqrt, linspace
from symcalc import construct_pdfs
from scipy.stats import lognorm
from sympy import lambdify
from scipy.integrate import trapz


def log_normal(x, d32, d03):
    sigma = log(d32 / d03)
    return lognorm.pdf(x, s=sigma, scale=d03)


class test_case:
    def __init__(
        self, nr_classes
    ):
        self.name = "c01"
        self.nr_classes = nr_classes
        dmean = 5.0e-4  # 500 microns
        dstd = 5.0e-5  # std of 10% from the mean
        dvar = dstd**2
        mu = log(dmean / sqrt(1 + dvar / (dmean**2)))
        sigma = sqrt(log(1 + dvar / dmean**2))
        v, volume_pdf, _, _ = construct_pdfs(mu, sigma)
        f = lambdify(v, volume_pdf, "numpy")
        vmax = 3e-10
        self.dv = vmax / nr_classes
        v = linspace(self.dv, vmax, nr_classes)
        print(trapz(f(v), x=v))
        print(trapz(v * f(v), x=v), dmean**3 * pi / 6)
        vmean = trapz(v * f(v), x=v)
        # Calculate total number of drops
        dpipe = 0.063
        NT = (0.062 * pi / 4 * dpipe**2) / vmean
        NTinlet = NT * 4.5 / 250
        # NTinlet = 0.062 / vmean
        self.Ninit = f(v) * self.dv * NTinlet


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

pbe_grids = [3, 10, 20, 50]
cases = [test_case(g) for g in pbe_grids]

if __name__ == "__main__":

    for c in cases:
        case_setup(c)
