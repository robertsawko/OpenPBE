#include <boost/python.hpp>
#include "../breakupKernels/CoulaloglouTavlarides.H"
#include "../coalescenceKernels/CoulaloglouTavlarides.H"

BOOST_PYTHON_MODULE(kernels)
{
    using namespace boost::python;
    using Foam::scalar;

    class_<Foam::breakupKernels::CoulaloglouTavlaridesImp>(
                "CoulaloglouTavlaridesBreakup",
                init<scalar, scalar, scalar, scalar>())
            .def("S", &Foam::breakupKernels::CoulaloglouTavlaridesImp::S);

    class_<Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl>(
                "CoulaloglouTavlaridesCoalescence",
                init<scalar, scalar, scalar, scalar>())
            .def("S", &Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl::S);
}
