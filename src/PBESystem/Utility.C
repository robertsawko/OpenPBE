#include "Utility.H"


void printAvgMaxMin(const Foam::fvMesh &mesh, const Foam::volScalarField &v)
{
    Foam::Info<< v.name() << ": avg, max,min "
        << v.weightedAverage(mesh.V()).value()
        << ", " << max(v).value()
        << ", " << min(v).value() << Foam::endl;
}
