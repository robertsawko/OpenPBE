from os import path
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
#import ipdb

templateCase = SolutionDirectory(
    "parameterized", archive=None, paraviewLink=False)

nr_classes = [10, 20, 40]

for nC in nr_classes:
    case = templateCase.cloneCase("testCase" + str(nC))
    phaseProperties = ParsedParameterFile(
        path.join(case.name, "constant", "phaseProperties"))
    phaseProperties["air"]["PBEDiameterCoeffs"]["MOCCoeffs"]["numberOfClasses"] = nC

    # manually fix bad pyfoam parsing
    phaseProperties["blending"]["default"]["type"] = "none"
    phaseProperties["drag"][1]["swarmCorrection"]["type"] = "none"

    phaseProperties.writeFile()
    n0 = ParsedParameterFile(path.join(case.name, "0", "n0"))
    for i in range(1, nC - 1):
        n0.header["object"] = "n" + str(i)
        n0.writeFileAs(path.join(case.name, "0", "n" + str(i)))
    n0["internalField"].setUniform(1)
    n0.header["object"] = "n" + str(nC - 1)
    n0.writeFileAs(path.join(case.name, "0", "n" + str(nC - 1)))

    #ipdb.set_trace()
    controlDict = ParsedParameterFile(
        path.join(case.name, "system", "controlDict")
    )
    controlDict["functions"]["probes"]["fields"] = [
        "n{0}".format(n) for n in range(nC)]
    controlDict.writeFile()
