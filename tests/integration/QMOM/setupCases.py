from os import path
from numpy import ones, array, arange
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


class test_case:
    def __init__(
        self, case_name, quadrature_order,
        initial_moments, phase_properties_name,
        end_time, delta_t
    ):
        self.name = case_name
        self.quadrature_order = quadrature_order
        self.number_of_moments = 2 * quadrature_order
        self.initial_moments = initial_moments
        self.phase_properties_name = phase_properties_name
        self.end_time = end_time
        self.delta_t = delta_t


class breakup_case(test_case):

    def __init__(self, quadrature_order):
        self.l = 1.0
        test_case.__init__(
            self, "pure_br", quadrature_order,
            initial_moments=ones(2 * quadrature_order),
            phase_properties_name="phaseProperties.breakup",
            end_time=10, delta_t=0.005
        )


class coalescence_case(test_case):

    def __init__(self, quadrature_order, N0=2, v0=0.5, C=0.1):
        self.N0 = N0
        self.v0 = v0
        self.C = C

        from sympy import simplify, symbols, integrate, exp, oo
        v = symbols('v', real=True, positive=True)

        initial_moments = 0.25 * N0 / v0**2 * array([
            float(simplify(integrate(v**(k + 1) * exp(-v / (2 * v0)), (v, 0, oo))))
            for k in range(2 * quadrature_order)])

        test_case.__init__(
            self, "pure_coal", quadrature_order,
            initial_moments,
            "phaseProperties.coalescence",
            1, 0.001
        )


def case_setup(ci):
    template_case = SolutionDirectory(
        "template", archive=None, paraviewLink=False)
    case = template_case.cloneCase(
        "{0}{1}".format(ci.name, ci.quadrature_order)
    )

    phase_properties = ParsedParameterFile(
        path.join("./diffs", ci.phase_properties_name))
    phase_properties["air"]["PBEDiameterCoeffs"]["QMOMCoeffs"]["quadratureOrder"] = ci.quadrature_order

    # manually fix bad pyfoam parsing
    phase_properties["blending"]["default"]["type"] = "none"
    phase_properties["drag"][1]["swarmCorrection"]["type"] = "none"
    phase_properties.writeFileAs(path.join(
        case.name, "constant", "phaseProperties"
    ))

    m0 = ParsedParameterFile(path.join(template_case.name, "0", "m0"))
    for i in range(ci.number_of_moments):
        m0.header["object"] = "m" + str(i)
        m0["internalField"].setUniform(ci.initial_moments[i])
        m0["dimensions"] = "[0 {0} 0 0 0 0 0]".format(3 * i)
        m0.writeFileAs(path.join(case.name, "0", "m" + str(i)))

    controlDict = ParsedParameterFile(
        path.join(case.name, "system", "controlDict")
    )
    controlDict["functions"]["probes"]["fields"] = [
        "m{0}".format(m) for m in range(ci.number_of_moments)]
    controlDict["endTime"] = ci.end_time
    controlDict["deltaT"] = ci.delta_t
    controlDict.writeFile()

pbe_grids = range(1, 5)
coalescence_cases = [coalescence_case(nC) for nC in pbe_grids]
breakup_cases = [breakup_case(nC) for nC in pbe_grids]

if __name__ == "__main__":

    for bc in breakup_cases:
        case_setup(bc)

    for cc in coalescence_cases:
        case_setup(cc)
